"""Cross-browser / cross-display legibility audit for the PeptiLine web app.

Loads every page under a matrix of display profiles (laptop, 27" monitor, tablet,
phone; device-pixel-ratios 1-3) in Chromium, Firefox and WebKit, and for each
page records:

  * the computed CSS font size of every visible piece of UI text, weighted by
    character count (min / p10 / median, and % of characters below 12 px and 14 px);
  * the same for text inside generated Plotly figures (SVG <text>);
  * horizontal page overflow (document wider than the viewport);
  * a full-page screenshot.

Not a pytest module (needs a running server and real browsers). Usage:

    bash dev.sh 8765 &
    python3 tests/display/display_audit.py --base-url http://127.0.0.1:8765 \
        --out /tmp/display_audit [--browsers chromium,firefox,webkit] \
        [--profiles laptop-15-scaled,desktop-27-qhd] [--generate]

Requires `pip install playwright && python3 -m playwright install` (WebKit on
Linux additionally needs `sudo python3 -m playwright install-deps webkit`).
"""

import argparse
import json
import statistics
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

# name: (css width, css height, device pixel ratio, is_mobile, what it stands for)
PROFILES = {
    "laptop-15-scaled": (1280, 720, 1.5, False, '15" 1080p laptop, 150% OS scaling'),
    "laptop-15": (1536, 864, 1.25, False, '15" 1080p laptop, 125% OS scaling'),
    "laptop-13-mac": (1440, 900, 2, False, '13-15" MacBook (Retina)'),
    "desktop-1080": (1920, 1080, 1, False, '24" 1080p monitor, 100%'),
    "desktop-27-qhd": (2560, 1440, 1, False, '27" QHD monitor, 100%'),
    "desktop-27-4k": (2560, 1440, 1.5, False, '27" 4K monitor, 150% OS scaling'),
    "tablet": (820, 1180, 2, True, "iPad, portrait"),
    "phone": (390, 844, 3, True, "iPhone-class phone"),
}

PAGES = {
    "landing": "/",
    "data_transformation": "/data_transformation/",
    "data_analysis": "/data_analysis/",
    "heatmap": "/heatmap/",
}

# Walks the DOM and returns [font_px, n_chars, tag, snippet] for each visible
# element that directly owns text. Plotly SVG text is reported separately.
MEASURE_JS = r"""
() => {
  const ui = [], plot = [];
  const visible = el => {
    const cs = getComputedStyle(el);
    if (cs.visibility === 'hidden' || cs.display === 'none' || +cs.opacity === 0) return false;
    const r = el.getBoundingClientRect();
    return r.width > 0 && r.height > 0;
  };
  const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  const seen = new Map();
  let n;
  while ((n = walker.nextNode())) {
    const t = n.textContent.replace(/\s+/g, ' ').trim();
    if (!t) continue;
    const el = n.parentElement;
    if (!el || ['SCRIPT', 'STYLE', 'NOSCRIPT', 'OPTION', 'TITLE'].includes(el.tagName)) continue;
    if (!visible(el)) continue;
    // skip text clipped away inside a closed/hidden ancestor (offsetParent null for HTML)
    if (!(el instanceof SVGElement) && el.offsetParent === null && getComputedStyle(el).position !== 'fixed') continue;
    seen.set(el, (seen.get(el) || '') + t);
  }
  for (const [el, t] of seen) {
    const px = parseFloat(getComputedStyle(el).fontSize);
    const rec = [px, t.length, el.tagName.toLowerCase(), t.slice(0, 50)];
    (el.closest('.js-plotly-plot') ? plot : ui).push(rec);
  }
  // <select> option text renders at the select's own size
  document.querySelectorAll('select').forEach(s => {
    if (!visible(s)) return;
    const px = parseFloat(getComputedStyle(s).fontSize);
    const txt = Array.from(s.options).slice(0, 5).map(o => o.text).join(' ');
    if (txt.trim()) ui.push([px, Math.min(txt.length, 80), 'select', txt.slice(0, 50)]);
  });
  const de = document.documentElement;
  const overflowX = Math.max(de.scrollWidth, document.body.scrollWidth) - window.innerWidth;
  const offenders = [];
  if (overflowX > 1) {
    document.querySelectorAll('body *').forEach(el => {
      const r = el.getBoundingClientRect();
      if (r.right > window.innerWidth + 1 && r.width > 0 && visible(el) && offenders.length < 8)
        offenders.push(el.tagName.toLowerCase() + (el.id ? '#' + el.id : '') +
          (el.className && typeof el.className === 'string' ? '.' + el.className.split(' ')[0] : '') +
          ' right=' + Math.round(r.right));
    });
  }
  return {ui, plot, overflowX, offenders,
          rootPx: parseFloat(getComputedStyle(de).fontSize),
          bodyPx: parseFloat(getComputedStyle(document.body).fontSize)};
}
"""


def summarize(records):
    if not records:
        return None
    sizes = []
    for px, n, *_ in records:
        sizes.extend([px] * n)
    sizes.sort()
    total = len(sizes)
    smallest = sorted(records, key=lambda r: r[0])[:5]
    return {
        "chars": total,
        "min_px": round(sizes[0], 1),
        "p10_px": round(sizes[int(total * 0.10)], 1),
        "median_px": round(statistics.median(sizes), 1),
        "pct_lt12": round(100 * sum(s < 12 for s in sizes) / total, 1),
        "pct_lt14": round(100 * sum(s < 14 for s in sizes) / total, 1),
        "smallest": [[round(r[0], 1), r[2], r[3]] for r in smallest],
    }


def wait_quiet(page, ms=600):
    try:
        page.wait_for_load_state("networkidle", timeout=15000)
    except Exception:
        pass
    page.wait_for_timeout(ms)


def populate(page, name, generate):
    """Load the bundled example data so selectors/settings panels are on screen."""
    if name == "data_transformation":
        for link in page.locator(".dt-load-example:visible").all()[:2]:
            link.click()
        wait_quiet(page, 1500)
    elif name == "data_analysis":
        page.click("#da-load-example")
        page.wait_for_function("!document.getElementById('generate-btn').disabled", timeout=30000)
        wait_quiet(page)
        if generate:
            page.click("#generate-btn")
            page.wait_for_selector(".js-plotly-plot, #plot-output img", timeout=120000)
            wait_quiet(page, 1500)
    elif name == "heatmap":
        page.click("#hm-load-example-csv")
        page.wait_for_timeout(1500)
        page.click("#hm-load-example-fasta")
        page.wait_for_function("!document.getElementById('generate-btn').disabled", timeout=30000)
        wait_quiet(page)
        if generate:
            # first protein + first two variables, as a user would pick them
            page.evaluate("""() => {
                for (const [id, k] of [['protein-selector', 1], ['varkey-selector', 2]]) {
                    const s = document.getElementById(id);
                    Array.from(s.options).slice(0, k).forEach(o => o.selected = true);
                    s.dispatchEvent(new Event('change', {bubbles: true}));
                }
            }""")
            page.wait_for_timeout(500)
            page.click("#generate-btn")
            page.wait_for_selector(".js-plotly-plot", state="visible", timeout=180000)
            wait_quiet(page, 2000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-url", default="http://127.0.0.1:8765")
    ap.add_argument("--out", default="display_audit_out")
    ap.add_argument("--browsers", default="chromium,firefox,webkit")
    ap.add_argument("--profiles", default=",".join(PROFILES))
    ap.add_argument("--pages", default=",".join(PAGES))
    ap.add_argument("--generate", action="store_true", help="also generate a plot on DA / heatmap")
    ap.add_argument("--zoom", type=float, default=1.0,
                    help="simulate browser zoom (e.g. 1.25) by shrinking the CSS viewport and raising DPR")
    args = ap.parse_args()

    out = Path(args.out)
    (out / "shots").mkdir(parents=True, exist_ok=True)
    results = []
    with sync_playwright() as p:
        for bname in args.browsers.split(","):
            try:
                browser = getattr(p, bname).launch()
            except Exception as e:
                print(f"!! {bname}: cannot launch ({str(e).splitlines()[0]})", file=sys.stderr)
                results.append({"browser": bname, "error": "launch failed"})
                continue
            for prof in args.profiles.split(","):
                w, h, dpr, mobile, desc = PROFILES[prof]
                z = args.zoom
                ctx_kwargs = dict(viewport={"width": int(w / z), "height": int(h / z)},
                                  device_scale_factor=dpr * z, has_touch=mobile)
                if bname != "firefox":  # Firefox does not support isMobile
                    ctx_kwargs["is_mobile"] = mobile
                ctx = browser.new_context(**ctx_kwargs)
                for pname in args.pages.split(","):
                    page = ctx.new_page()
                    rec = {"browser": bname, "profile": prof, "desc": desc, "page": pname, "zoom": z}
                    try:
                        page.goto(args.base_url + PAGES[pname], wait_until="load", timeout=60000)
                        wait_quiet(page)
                        if pname != "landing":
                            populate(page, pname, args.generate)
                        m = page.evaluate(MEASURE_JS)
                        rec.update(ui=summarize(m["ui"]), plot=summarize(m["plot"]),
                                   overflow_x=m["overflowX"], overflow_offenders=m["offenders"],
                                   root_px=m["rootPx"], body_px=m["bodyPx"])
                        shot = out / "shots" / f"{bname}_{prof}_z{z}_{pname}.png"
                        page.screenshot(path=str(shot), full_page=True)
                        rec["shot"] = str(shot)
                    except Exception as e:
                        rec["error"] = str(e).splitlines()[0][:200]
                    ui = rec.get("ui") or {}
                    print(f"{bname:9s} {prof:17s} {pname:20s} root={rec.get('root_px')} "
                          f"min={ui.get('min_px')} p10={ui.get('p10_px')} med={ui.get('median_px')} "
                          f"<12px={ui.get('pct_lt12')}% overflowX={rec.get('overflow_x')} "
                          f"{rec.get('error', '')}", flush=True)
                    results.append(rec)
                    page.close()
                ctx.close()
            browser.close()
    (out / "results.json").write_text(json.dumps(results, indent=1))
    print(f"\nwrote {out / 'results.json'}")


if __name__ == "__main__":
    main()
