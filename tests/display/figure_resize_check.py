"""Check that generated figures re-flow across display sizes.

Generates the default example figure on the Data Analysis and Heat Map pages,
then resizes the browser window through laptop / 27" monitor / tablet / phone
widths. At each size it records the rendered Plotly width vs. its container,
the smallest on-screen plot text, and saves a screenshot of the figure.

    python3 tests/display/figure_resize_check.py --base-url http://127.0.0.1:8765 --out /tmp/fig_check
"""

import argparse
import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

sys.path.insert(0, str(Path(__file__).parent))
from display_audit import populate, wait_quiet  # noqa: E402

SIZES = [(1280, 720), (1920, 1080), (2560, 1440), (820, 1180), (390, 844)]

MEASURE = r"""
() => Array.from(document.querySelectorAll('.js-plotly-plot'))
  .filter(p => p.getBoundingClientRect().height > 0)
  .map(p => {
    const svg = p.querySelector('.main-svg');
    const texts = Array.from(p.querySelectorAll('text'))
      .filter(t => t.textContent.trim() && t.getBoundingClientRect().width > 0)
      .map(t => parseFloat(getComputedStyle(t).fontSize));
    return {id: p.id, container_w: Math.round(p.parentElement.getBoundingClientRect().width),
            plot_w: Math.round(svg ? svg.getBoundingClientRect().width : 0),
            min_text_px: texts.length ? Math.min(...texts) : null,
            n_text: texts.length,
            page_overflow: document.documentElement.scrollWidth - innerWidth};
  })
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-url", default="http://127.0.0.1:8765")
    ap.add_argument("--out", default="figure_check_out")
    ap.add_argument("--browsers", default="chromium,firefox,webkit")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    results = []
    with sync_playwright() as p:
        for bname in args.browsers.split(","):
            try:
                browser = getattr(p, bname).launch()
            except Exception as e:
                print(f"!! {bname}: cannot launch ({str(e).splitlines()[0]})")
                continue
            for pname, path in (("data_analysis", "/data_analysis/"), ("heatmap", "/heatmap/")):
                page = browser.new_page(viewport={"width": SIZES[0][0], "height": SIZES[0][1]})
                page.goto(args.base_url + path)
                wait_quiet(page)
                populate(page, pname, generate=True)
                for w, h in SIZES:
                    page.set_viewport_size({"width": w, "height": h})
                    page.wait_for_timeout(1200)  # Plotly's responsive resize is debounced
                    figs = page.evaluate(MEASURE)
                    for f in figs:
                        fits = f["plot_w"] <= f["container_w"] + 2
                        # a wide figure is fine if it scrolls inside its own frame
                        status = "OK" if fits else ("SCROLLS-IN-FRAME" if f["page_overflow"] <= 1 else "PAGE-OVERFLOW")
                        print(f"{bname:9s} {pname:14s} {w:5d}px  {f['id'] or '-':28s} plot={f['plot_w']}/"
                              f"{f['container_w']} {status} min_text={f['min_text_px']}px")
                        results.append(dict(browser=bname, page=pname, width=w, **f, status=status))
                    loc = page.locator(".js-plotly-plot >> visible=true").first
                    if loc.count():
                        loc.screenshot(path=str(out / f"{bname}_{pname}_{w}.png"))
                page.close()
            browser.close()
    (out / "results.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
