"""Draw the phase stack and two existing TESS deficit images."""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def draw(audit, pixels):
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 960 490" role="img">',
           '<title>TIC 238187210 transit stack and TESS deficit maps</title>',
           '<rect width="960" height="490" fill="#f7f9fc"/>',
           '<g font-family="system-ui, sans-serif" fill="#182937">']

    def label(x, y, words, size=14, anchor="start"):
        svg.append(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}">{words}</text>')

    label(36, 42, "TIC 238187210", 26)
    label(36, 67, "TESS light curves and in-transit deficit images", 16)
    label(40, 113, "Phase stack at 26.02685 days", 16)
    stacks = audit["stacks"]
    hour = np.concatenate([event["hour"] for event in stacks])
    flux = 100 * np.concatenate([event["relative"] for event in stacks])

    def xy(t, f):
        return 87 + (t + 8) * 26, 139 + (.5 - f) * 125

    for f in (.5, 0, -.5, -1, -1.5, -2):
        x, y = xy(-8, f)
        svg.append(f'<path d="M{x} {y}h416" stroke="#ccd5dd"/>')
        label(x - 9, y + 5, f"{f:g}%", 13, "end")
    for t in (-8, -4, 0, 4, 8):
        x, y = xy(t, -2)
        label(x, y + 20, str(t), 13, "middle")
    for t, f in zip(hour, flux):
        if abs(t) <= 8 and -2 <= f <= .5:
            x, y = xy(t, f)
            svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1" fill="#65859b" opacity=".18"/>')
    for lo in np.arange(-8, 8, .4):
        values = flux[(hour >= lo) & (hour < lo + .4)]
        if len(values) >= 5:
            x, y = xy(lo + .2, np.median(values))
            svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="#075570"/>')
    label(295, 447, "Hours from predicted midpoint", 14, "middle")
    for i, panel in enumerate(pixels):
        left = 550 + 195 * i
        top = 139
        tile = 25
        label(left, 113, f'Sector {panel["sector"]}', 16)
        values = np.array(panel["missing"])
        for (row, col), val in np.ndenumerate(values):
            tint = min(1, max(0, val / 3))
            red, green, blue = (round(a + tint * (b - a)) for a, b in zip((231, 238, 245), (17, 71, 95)))
            svg.append(f'<rect x="{left + col * tile}" y="{top + (6-row) * tile}" width="{tile}" height="{tile}" fill="rgb({red},{green},{blue})" stroke="#f7f9fc"/>')
        stars = {str(star["tic"]): star for star in panel["catalog_pixels"]}
        for tic, color in (("238187210", "#f1a52c"), ("238187209", "#d02c7a")):
            star = stars[tic]
            x = left + (star["x"] + .5) * tile
            y = top + (6.5 - star["y"]) * tile
            svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="7" fill="none" stroke="{color}" stroke-width="3"/>')
    label(550, 343, "Darker pixels show larger flux deficits.", 13)
    for y, color, words in ((368, "#f1a52c", "Target TIC 238187210"),
                            (391, "#d02c7a", "Bright neighbor TIC 238187209")):
        svg.append(f'<circle cx="558" cy="{y - 4}" r="6" fill="none" stroke="{color}" stroke-width="3"/>')
        label(572, y, words, 13)
    label(36, 474, "Dots show 24-minute bin medians. The stack and deficit maps share TESS observations.", 13)
    return "\n".join(svg + ["</g></svg>"])


if __name__ == "__main__":
    import sys

    audit = json.loads((HERE / "tic238_audit.json").read_text(encoding="utf-8"))
    pixels = json.loads((HERE / "tic238_pixel_inputs.json").read_text(encoding="utf-8"))
    Path(sys.argv[1]).write_text(draw(audit, pixels), encoding="utf-8")
