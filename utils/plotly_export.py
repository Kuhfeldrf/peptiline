"""Server-side static-image export for Plotly figures via Kaleido.

Renders publication-grade files from a figure's JSON on the server, so the
downloaded output is resolution-deterministic — independent of the viewer's
browser window (unlike client-side ``Plotly.downloadImage``, whose pixel size
tracks the on-screen element). Shared by the Data Analysis and Heatmap
dashboards.

Formats:
  * PNG — rasterized at ``PUBLICATION_DPI`` with the DPI written into the file's
    metadata (via Pillow), so a journal's automated resolution check reads it.
  * SVG — true vector output; scales losslessly, text stays text.

Kaleido v1+ downloads its own headless Chrome build (see Dockerfile's
`kaleido_get_chrome` step) rather than bundling one, so no system Chrome
package is needed, but the Chrome binary's shared-library dependencies
(libnss3, libgbm1, etc. — also installed in the Dockerfile) must be present.

Browser reuse: left to itself, ``fig.to_image()`` launches a fresh headless
Chrome, renders one figure, and tears it back down -- ~3s of pure process
overhead per export, on top of the render itself. ``_ensure_kaleido_server``
starts Kaleido's persistent browser server on first use and leaves it running
for the life of the process, so every export after the first in a given
gunicorn worker skips the launch/teardown entirely. Safe to call repeatedly
and from multiple threads (gunicorn runs this app with gthread workers).
"""
import atexit
import io
import threading
import warnings

import kaleido
import plotly.io as pio
from PIL import Image

# plotly.io.to_image always builds a non-empty `kopts` dict (its default
# `headers` is `{'X-Requested-With': 'plotly.py'}`, never empty), and
# kaleido's persistent-server path warns on every single call when kopts is
# non-empty -- so once the server is running, this fires on every export.
# Nothing we pass through kopts (we don't set plotlyjs/mathjax/headers
# overrides) actually depends on per-call kopts once the server is up; it's
# only meaningful at server-start time. Scoped to this exact message so it
# doesn't mask an unrelated warning.
warnings.filterwarnings(
    'ignore',
    message='The kopts argument is ignored if using a server.',
    category=UserWarning,
)

_kaleido_server_lock = threading.Lock()
_kaleido_server_started = False


def _ensure_kaleido_server():
    global _kaleido_server_started
    if _kaleido_server_started:
        return
    with _kaleido_server_lock:
        if _kaleido_server_started:
            return
        try:
            kaleido.start_sync_server(silence_warnings=True)
            atexit.register(kaleido.stop_sync_server, silence_warnings=True)
        except Exception:
            # Fall back to kaleido's default one-shot-per-call behavior --
            # slower (pays the browser launch cost every export) but still
            # correct. Nothing here is essential to the export succeeding.
            pass
        _kaleido_server_started = True

# Pillow's decompression-bomb heuristic assumes an untrusted/unknown source;
# these PNGs are rendered server-side from our own figures, and a large
# publication figure at PUBLICATION_DPI routinely exceeds Pillow's default
# 89-megapixel ceiling. Raise it well above what any real figure needs
# instead of disabling the check outright.
Image.MAX_IMAGE_PIXELS = 300_000_000

# Publication resolution. The target journal (JPR/ACS) specifies a minimum of
# 300 dpi for COLOUR art (docs/JPR_journal_requirements.md); 600 dpi is used
# here to sit comfortably above that floor and also clear the 600-dpi grayscale
# threshold, giving headroom for any figure type.
PUBLICATION_DPI = 600

# Plotly lays figures out in CSS pixels at a nominal 96 px/inch; scaling the
# raster by dpi/96 yields the requested print DPI.
_CSS_PPI = 96


def _figure(figure_json):
    """Build a Plotly figure from a JSON string (the transport both dashboards
    already hold client-side)."""
    if not isinstance(figure_json, str):
        # Accept a dict too, for callers that pass the parsed figure.
        import json
        figure_json = json.dumps(figure_json)
    return pio.from_json(figure_json)


def png_bytes(figure_json, width=None, height=None, dpi=PUBLICATION_DPI):
    """Return PNG bytes at ``dpi``, with the resolution stamped in metadata.

    ``width``/``height`` (CSS px, e.g. the on-screen figure size) fix the print
    dimensions; omit them to use the figure's own layout size.
    """
    _ensure_kaleido_server()
    fig = _figure(figure_json)
    scale = dpi / _CSS_PPI
    raw = fig.to_image(format='png', width=width, height=height, scale=scale)
    # Re-save through Pillow purely to embed the DPI (pHYs) metadata; pixels
    # are unchanged.
    im = Image.open(io.BytesIO(raw))
    out = io.BytesIO()
    im.save(out, format='PNG', dpi=(dpi, dpi))
    return out.getvalue()


def svg_bytes(figure_json, width=None, height=None):
    """Return true-vector SVG bytes (resolution-independent)."""
    _ensure_kaleido_server()
    fig = _figure(figure_json)
    return fig.to_image(format='svg', width=width, height=height)
