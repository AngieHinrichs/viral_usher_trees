#!/usr/bin/env python3

import sys
import io
import html as _html
import xml.etree.ElementTree as ET

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------
# Input files
# ----------------------------

TREE_FILE = "../tree_metadata.tsv"
BALT_FILE = "baltimore_classification.tsv"
OUTPUT = "baltimore_classification_plot.pdf"
HTML_OUTPUT = "baltimore_classification_plot.html"

HTML_MODE = True

# ----------------------------
# Read data
# ----------------------------

trees = pd.read_csv(TREE_FILE, sep="\t")
balt = pd.read_csv(BALT_FILE, sep="\t")

# Merge
df = trees.merge(
    balt[["Sample", "Baltimore_group"]],
    left_on="tree_name",
    right_on="Sample",
    how="left"
)

# Desired order
order = ["I", "II", "III", "IV", "V", "VI", "VII", "viroid"]

missing_count = df["Baltimore_group"].isna().sum()
if missing_count:
    df["Baltimore_group"] = df["Baltimore_group"].fillna("Unclassified")
    order.append("Unclassified")

df["Baltimore_group"] = pd.Categorical(
    df["Baltimore_group"],
    categories=order,
    ordered=True
)

# ----------------------------
# Summary statistics
# ----------------------------

summary = (
    df.groupby("Baltimore_group", observed=False)
      .agg(
          total_tips=("tip_count", "sum"),
          n_trees=("tip_count", "size")
      )
      .reindex(order)
      .fillna(0)
)

# ----------------------------
# Plot
# ----------------------------

fig = plt.figure(figsize=(12,6))

gs = fig.add_gridspec(1,2, width_ratios=[1,2])

###########################################################
# Panel A
###########################################################

ax1 = fig.add_subplot(gs[0])

bars = ax1.bar(
    summary.index,
    summary["total_tips"],
)

for _j, _bar in enumerate(bars):
    _bar.set_gid(f"barA_{_j}")

ax1.set_yscale("log")
ax1.set_ylabel("Total genomes (tips)")
ax1.set_xlabel("Baltimore class")
ax1.set_title("A")

###########################################################
# Panel B
###########################################################

ax2 = fig.add_subplot(gs[1])

x_positions = np.arange(len(order))

rng = np.random.default_rng(1)

_hover_rows = {}

for i, group in enumerate(order):

    subset = df[df["Baltimore_group"] == group]

    if len(subset) == 0:
        continue

    jitter = rng.normal(0, 0.08, len(subset))

    ax2.scatter(
        np.full(len(subset), i) + jitter,
        subset["tip_count"],
        s=18,
        alpha=0.7,
    )

    ax2.collections[-1].set_gid(f"grpB_{i}")
    _hover_rows[i] = subset

ax2.set_yscale("log")
ax2.set_xticks(x_positions)
ax2.set_xticklabels(order)
ax2.set_ylabel("Tips per tree")
ax2.set_xlabel("Baltimore class")
ax2.set_title("B")

plt.tight_layout()

plt.savefig(OUTPUT)

# ----------------------------
# Optional interactive HTML (--html)
# ----------------------------

if HTML_MODE:

    SVGNS = "http://www.w3.org/2000/svg"
    ET.register_namespace("", SVGNS)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

    _buf = io.StringIO()
    plt.savefig(_buf, format="svg")
    _root = ET.fromstring(_buf.getvalue())
    _by_id = {e.get("id"): e for e in _root.iter() if e.get("id")}

    # fields shown on hover, in order; blanks are skipped per-row
    _fields = ["accession", "organism", "isolate", "strain",
               "serotype", "segment", "Taxonomy_ID", "tip_count"]

    def _tip(pairs):
        out = []
        for k, v in pairs:
            if v is None:
                continue
            v = str(v).strip()
            if v == "" or v.lower() in ("nan", "none"):
                continue
            out.append(f"<b>{_html.escape(str(k))}</b> {_html.escape(v)}")
        return "<br>".join(out)

    def _mark(el, tiphtml):
        if el is None or not tiphtml:
            return
        el.set("class", "hot")
        el.set("data-tip", tiphtml)

    # --- Panel A: one tooltip per bar ---
    for _j, _grp in enumerate(order):
        _row = summary.loc[_grp]
        _mark(_by_id.get(f"barA_{_j}"), _tip([
            ("Baltimore class", _grp),
            ("Total genomes (tips)", f"{int(_row['total_tips']):,}"),
            ("Trees", int(_row["n_trees"])),
        ]))

    # --- Panel B: one tooltip per point ---
    # <use> elements inside a PathCollection <g> follow data order
    for _i, _subset in _hover_rows.items():
        _g = _by_id.get(f"grpB_{_i}")
        if _g is None:
            continue
        _uses = list(_g.iter(f"{{{SVGNS}}}use"))
        for _use, (_, _r) in zip(_uses, _subset.iterrows()):
            _pairs = [("tree", _r.get("tree_name")),
                      ("Baltimore class", order[_i])]
            _pairs += [(_f, _r.get(_f)) for _f in _fields]
            _mark(_use, _tip(_pairs))

    _svg = ET.tostring(_root, encoding="unicode")

    _page = """<!DOCTYPE html>
<meta charset="utf-8">
<title>resource_summary</title>
<style>
  body { margin: 0; padding: 24px; font: 13px/1.45 -apple-system, Segoe UI, Helvetica, Arial, sans-serif; }
  #wrap svg { max-width: 100%; height: auto; }
  .hot { cursor: pointer; }
  .hot:hover { outline: none; filter: brightness(0.55); }
  #tip {
    position: fixed; z-index: 10; display: none; pointer-events: none;
    max-width: 340px; padding: 7px 10px;
    background: rgba(20,20,20,.94); color: #fff; border-radius: 5px;
    font-size: 12px; line-height: 1.5; box-shadow: 0 2px 10px rgba(0,0,0,.3);
  }
  #tip b { color: #9ecbff; font-weight: 600; }
</style>
<div id="wrap">__SVG__</div>
<div id="tip"></div>
<script>
(function () {
  var tip = document.getElementById('tip');
  var wrap = document.getElementById('wrap');

  function place(e) {
    var pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
    var x = e.clientX + pad, y = e.clientY + pad;
    if (x + w > window.innerWidth)  x = e.clientX - w - pad;
    if (y + h > window.innerHeight) y = e.clientY - h - pad;
    tip.style.left = Math.max(0, x) + 'px';
    tip.style.top  = Math.max(0, y) + 'px';
  }

  wrap.addEventListener('pointermove', function (e) {
    var el = e.target.closest('.hot');
    if (!el) { tip.style.display = 'none'; return; }
    tip.innerHTML = el.getAttribute('data-tip');
    tip.style.display = 'block';
    place(e);
  });

  wrap.addEventListener('pointerleave', function () {
    tip.style.display = 'none';
  });
})();
</script>
"""
    with open(HTML_OUTPUT, "w") as _fh:
        _fh.write(_page.replace("__SVG__", _svg))

plt.close()

print(summary)
print(f"\nSaved {OUTPUT}")
if HTML_MODE:
    print(f"Saved {HTML_OUTPUT}")
