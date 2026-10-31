from pathlib import Path
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "figures"
OUT.mkdir(exist_ok=True)

from figstyle import *
import pandas as pd, numpy as np
R = str(REPO / "results" / "study2" / "main") + "/"
LAB = {"Cox+DCA − Clinical-Cox": "Cox+DCA vs clinical Cox (primary)",
       "Cox+DCA − Cox+Attention": "Cox+DCA vs Cox+attention",
       "Cox+DCA − Cox+Concat": "Cox+DCA vs Cox+concatenation",
       "Cox+DCA − Cox+Gated": "Cox+DCA vs Cox+gated",
       "Cox+DCA − Stacked-Cox": "Cox+DCA vs stacked Cox",
       "Cox+DCA − Late-Fusion": "Cox+DCA vs late fusion",
       "Cox+WSI − Clinical-Cox": "Cox+WSI vs clinical Cox",
       "Cox+Omics − Clinical-Cox": "Cox+omics vs clinical Cox",
       "DCA − Attention-Fusion": "DCA vs attention fusion",
       "DCA − Concat-Fusion": "DCA vs concatenation fusion",
       "DCA − Gated-Fusion": "DCA vs gated fusion"}
fig, axes = plt.subplots(1, 2, figsize=(173 * MM, 78 * MM), sharey=True)
fig.subplots_adjust(left=0.27, right=0.99, top=0.88, bottom=0.14, wspace=0.08)
for ax, c, t in zip(axes, ["STAD", "CRC"], ["A", "B"]):
    d = pd.read_csv(R + f"v2_{c}_uni_cap500_contrasts.csv")
    d["lab"] = d.contrast.map(LAB)
    for i, r in enumerate(d.itertuples()):
        col = RED if r.ci_hi < 0 else "black"
        ax.errorbar(r.diff_c, i, xerr=[[r.diff_c - r.ci_lo], [r.ci_hi - r.diff_c]], fmt="D" if i == 0 else "o",
                    color=col, ms=4.5 if i == 0 else 3.5, elinewidth=1.4 if i == 0 else 0.8, capsize=1.5)
        p = 2 * min(r.share_boot_le_0, 1 - r.share_boot_le_0)
        f = lambda v: ("0.000" if abs(v) < 0.0005 and v != r.ci_hi else ("−0.0002" if abs(v) < 0.0005 else f"{v:+.3f}".replace("-", "−")))
        ax.text(1.0, i, f"{f(r.diff_c)} ({f(r.ci_lo)}, {f(r.ci_hi)})", transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=7.2, color=col)
    ax.axvline(0, color=GRAY, lw=0.7, ls="--")
    ax.axhspan(-0.5, 0.5, color=BLUE, alpha=0.08, lw=0)
    ax.axhline(5.5, color=GRAY, lw=0.4); ax.axhline(7.5, color=GRAY, lw=0.4)
    ax.set_xlim(-0.08, 0.28); ax.set_xticks([-0.05, 0, 0.05, 0.10])
    ax.set_xlabel("Difference in pooled C-index (95% bootstrap CI)")
    ax.set_title({"STAD": "Gastric cancer (278 patients, 108 deaths)", "CRC": "Colorectal cancer (358 patients, 74 deaths)"}[c], fontsize=8)
    ax.text(-0.02 if c == "CRC" else -0.36, 1.07, t, transform=ax.transAxes, fontsize=8, fontweight="bold")
    ax.tick_params(axis="y", length=0)
axes[0].set_yticks(range(len(d))); axes[0].set_yticklabels(d.lab); axes[0].invert_yaxis()
axes[0].get_yticklabels()[0].set_fontweight("bold")
fig.savefig(str(OUT / "study2_contrasts.png"), dpi=600); fig.savefig(str(OUT / "study2_contrasts.tif"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
print("ok")
