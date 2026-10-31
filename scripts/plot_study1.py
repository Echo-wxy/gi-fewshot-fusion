from pathlib import Path
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "figures"
OUT.mkdir(exist_ok=True)

from figstyle import *
import pandas as pd, numpy as np
R = str(REPO / "results" / "study1") + "/"
s = pd.read_csv(R + "summary.csv"); p = pd.read_csv(R + "paired_differences.csv")
mp = pd.read_csv(R + "masked_paired_differences.csv")
NAMES = {"NC-L2N": "Nearest centroid", "SimpleShot-CL2N": "SimpleShot", "LogReg-CL2N": "Logistic regression",
         "NonEpisodic-Head+NC": "Non-episodic adapter", "ProtoNet-Head (ours)": "Episodic adapter (proposed)"}
STY = {"NC-L2N": dict(color=GRAY, ls=":", marker="v", lw=1.0, mfc="white"),
       "SimpleShot-CL2N": dict(color=GRAY, ls="--", marker="s", lw=1.0, mfc="white"),
       "LogReg-CL2N": dict(color=GRAY, ls="-.", marker="D", lw=1.0, mfc="white", ms=3.5),
       "NonEpisodic-Head+NC": dict(color=TEAL, ls="--", marker="o", lw=1.0, mfc="white", alpha=0.8),
       "ProtoNet-Head (ours)": dict(color=TEAL_E, ls="-", marker="o", lw=2.0)}
fig = plt.figure(figsize=(173 * MM, 66 * MM))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 0.9, 1.25], wspace=0.42, left=0.065, right=0.99, bottom=0.19, top=0.74)
# (a) accuracy vs shot
ax = fig.add_subplot(gs[0])
m = s[(s.split == "main_rare5") & (s.sampling == "scene")]
for k, st in STY.items():
    g = m[m.method == k].sort_values("shot")
    yerr = np.vstack([g.acc - g.acc_cluster_ci_lo, g.acc_cluster_ci_hi - g.acc])
    ax.errorbar(g.shot + {"NC-L2N": -0.3, "SimpleShot-CL2N": -0.15, "LogReg-CL2N": 0, "NonEpisodic-Head+NC": 0.15,
                          "ProtoNet-Head (ours)": 0.3}[k], g.acc, yerr=yerr, capsize=1.5, elinewidth=0.6,
                ms=st.pop("ms", 4), label=NAMES[k], **st)
ax.set_xticks([1, 5, 10]); ax.set_xlabel("Support images per class"); ax.set_ylabel("Accuracy (%)")
ax.set_ylim(35, 70)
# (b) paired gains, main split
ax2 = fig.add_subplot(gs[1])
x = np.arange(3); w = 0.2
for j, (vs, col, lab) in enumerate([("SimpleShot-CL2N", GRAY, "vs SimpleShot"), ("NonEpisodic-Head+NC", TEAL, "vs non-episodic adapter")]):
    for src, off, hatch, alpha in ((p, -w / 2, None, 1.0), (mp, w / 2, "////", 0.45)):
        g = src[(src.split == "main_rare5") & (src.sampling == "scene") & (src.versus == vs)].sort_values("shot")
        xx = x + (j - 0.5) * 2 * w + off
        cols = [RED if hi < 0 else col for hi in g.ci_hi]
        ax2.bar(xx, g.diff_pp, w, color=cols, alpha=alpha, hatch=hatch, edgecolor=cols, lw=0.5)
        ax2.errorbar(xx, g.diff_pp, yerr=[g.diff_pp - g.ci_lo, g.ci_hi - g.diff_pp], fmt="none", ecolor="black", elinewidth=0.6, capsize=1.2)
ax2.axhline(0, color="black", lw=0.6)
ax2.set_xticks(x); ax2.set_xticklabels(["1", "5", "10"]); ax2.set_xlabel("Support images per class")
ax2.set_ylabel("Paired gain (pp)")
# (c) 5-shot gain across splits
ax3 = fig.add_subplot(gs[2])
order = [("main_rare5", "scene", "Rare lesions (main)"), ("family_UC4", "scene", "UC grades"), ("family_esophagus4", "scene", "Esophageal lesions"),
         ("random5_1", "scene", "Random 1"), ("random5_2", "strict", "Random 2*"), ("random5_3", "scene", "Random 3"),
         ("random5_4", "scene", "Random 4"), ("random5_5", "scene", "Random 5")]
for i, (sp, md, lab) in enumerate(order):
    for vs, col, dy, mk in (("SimpleShot-CL2N", GRAY, -0.15, "s"), ("NonEpisodic-Head+NC", TEAL, 0.15, "o")):
        r = p[(p.split == sp) & (p.sampling == md) & (p.shot == 5) & (p.versus == vs)].iloc[0]
        c = RED if r.ci_hi < 0 else col
        ax3.errorbar(r.diff_pp, i + dy, xerr=[[r.diff_pp - r.ci_lo], [r.ci_hi - r.diff_pp]], fmt=mk, color=c, ms=3.5,
                     mfc=c, elinewidth=0.7, capsize=1.2)
ax3.axvline(0, color="black", lw=0.6)
ax3.set_yticks(range(len(order))); ax3.set_yticklabels([o[2] for o in order]); ax3.invert_yaxis()
ax3.set_xlabel("5-shot paired gain of proposed adapter (pp)")
ax3.tick_params(axis="y", length=0)
# legends on top
h, l = ax.get_legend_handles_labels()
fig.legend(h, l, loc="upper left", bbox_to_anchor=(0.0, 0.995), ncol=3, frameon=False, handlelength=2.2, columnspacing=1.0)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
fig.legend([Patch(fc=GRAY, ec=GRAY), Patch(fc=TEAL, ec=TEAL), Patch(fc="white", ec="black", hatch="////", lw=0.5)],
           ["vs SimpleShot", "vs non-episodic adapter", "overlay masked"], loc="upper right", bbox_to_anchor=(1.0, 0.995),
           ncol=1, frameon=False, handlelength=1.4, labelspacing=0.25)
for a, t in ((ax, "A"), (ax2, "B"), (ax3, "C")):
    a.text(-0.2 if a is not ax3 else -0.42, 1.04, t, transform=a.transAxes, fontsize=8, fontweight="bold", va="bottom")
fig.savefig(str(OUT / "study1_fewshot.png"), dpi=600); fig.savefig(str(OUT / "study1_fewshot.tif"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
print("ok")
