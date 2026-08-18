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

Kaleido 0.2.x ships its own headless Chromium, so this adds no separate
system-Chrome dependency on the server.
"""
import io

import plotly.io as pio
from PIL import Image

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
    fig = _figure(figure_json)
    return fig.to_image(format='svg', width=width, height=height)
