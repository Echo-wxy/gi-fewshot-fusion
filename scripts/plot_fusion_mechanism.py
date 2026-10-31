from pathlib import Path
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "figures"
OUT.mkdir(exist_ok=True)

from figstyle import *
import pandas as pd, numpy as np
R = str(REPO / "results" / "study2" / "main") + "/"
MOD = [("wsi", "WSI", PURPLE), ("omics", "Omics", AMBER), ("clin", "Clinical", BLUE)]
fig = plt.figure(figsize=(173 * MM, 64 * MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.35], wspace=0.28, left=0.07, right=0.99, bottom=0.2, top=0.8)
ax = fig.add_subplot(gs[0])
pos = 0; ticks = []; tl = []
for c, lab in (("STAD", "Gastric"), ("CRC", "Colorectal")):
    w = pd.read_csv(R + f"v2_{c}_uni_cap500_fusion_weights.csv"); w = w[w.method == "Cox+DCA"]
    for j, (k, name, col) in enumerate(MOD):
        bp = ax.boxplot(w[f"alpha_{k}"], positions=[pos + j], widths=0.6, patch_artist=True, showfliers=False,
                        medianprops=dict(color=col, lw=1.6), whiskerprops=dict(color=col, lw=0.8),
                        capprops=dict(color=col, lw=0.8), boxprops=dict(facecolor=col, alpha=0.25, edgecolor=col, lw=0.8))
    ticks.append(pos + 1); tl.append(lab); pos += 4
ax.axhline(1 / 3, color=GRAY, lw=0.6, ls="--")
ax.set_xticks(ticks); ax.set_xticklabels(tl); ax.set_ylabel("Fusion weight per patient (Cox+DCA)")
ax.set_ylim(0, 1); ax.tick_params(axis="x", length=0)
ax2 = fig.add_subplot(gs[1])
rng = np.random.default_rng(0)
pos = 0; ticks = []; tl = []
for c, lab in (("STAD", "Gastric"), ("CRC", "Colorectal")):
    m = pd.read_csv(R + f"v2_{c}_uni_cap500_fusion_mechanism.csv")
    for meth, mk in (("DCA", "o"), ("Cox+DCA", "D")):
        g = m[m.method == meth]
        for j, (k, name, col) in enumerate(MOD):
            v = g[f"dNLL_remove_{k}"].to_numpy(); x = pos + j * 0.28
            ax2.scatter(x + rng.uniform(-0.07, 0.07, len(v)), v, s=7, color=col, alpha=0.55, marker=mk, lw=0)
            ax2.plot([x - 0.11, x + 0.11], [np.median(v)] * 2, color=col, lw=1.8)
        ticks.append(pos + 0.28); tl.append(f"{lab}\n{meth}"); pos += 1.1
    pos += 0.3
ax2.axhline(0, color="black", lw=0.6)
ax2.set_xticks(ticks); ax2.set_xticklabels(tl); ax2.tick_params(axis="x", length=0)
ax2.set_ylabel("Change in test NLL when modality\nis replaced by its training mean")
ax2.set_ylim(-0.32, 0.30)
from matplotlib.patches import Patch
fig.legend([Patch(fc=c_, ec=c_, alpha=0.6) for _, _, c_ in MOD], [n for _, n, _ in MOD], loc="upper center",
           bbox_to_anchor=(0.5, 0.995), ncol=3, frameon=False)
ax.text(-0.2, 1.05, "A", transform=ax.transAxes, fontsize=8, fontweight="bold")
ax2.text(-0.17, 1.05, "B", transform=ax2.transAxes, fontsize=8, fontweight="bold")
fig.savefig(str(OUT / "study2_fusion_mechanism.png"), dpi=600); fig.savefig(str(OUT / "study2_fusion_mechanism.tif"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
print("ok")
