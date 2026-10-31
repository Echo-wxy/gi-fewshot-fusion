from pathlib import Path
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "figures"
OUT.mkdir(exist_ok=True)

from figstyle import *
import pandas as pd, numpy as np
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
P = str(REPO / "results" / "study2") + "/"
O = P + "main/"
res = []
fig, axes = plt.subplots(2, 2, figsize=(173 * MM, 146 * MM), sharey=True)
fig.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.15, hspace=0.78, wspace=0.26)
HI, LO = "#2B2B2B", "#9A9A9A"
ticks = [0, 12, 24, 36, 48, 60]
groups = {}
for r, c in enumerate(["STAD", "CRC"]):
    tab = pd.read_csv(P + f"cohort_{c}.csv").set_index("patient")
    oof = pd.read_csv(O + f"v2_{c}_uni_cap500_oof_predictions.csv")
    for k, (meth, lab) in enumerate([("Clinical-Cox", "Clinical Cox"), ("Cox+DCA", "Cox+DCA")]):
        rk = oof[oof.method == meth].groupby("patient").risk_rank_in_fold.mean()
        hi = (rk.reindex(tab.index) > 0.5).astype(int)
        groups[(c, meth)] = hi
        T, E = tab.time.values, tab.event.values
        ax = axes[r, k]
        for g, col, ls, name in ((0, LO, "--", "Low risk"), (1, HI, "-", "High risk")):
            m = hi.values == g
            kmf = KaplanMeierFitter().fit(T[m], E[m], label=name)
            kmf.plot_survival_function(ax=ax, ci_show=False, color=col, ls=ls, lw=1.4,
                                       show_censors=True, censor_styles={"ms": 3, "marker": "|"})
            # number at risk
            for t in ticks:
                ax.text(t, -0.36 - 0.09 * g, f"{int((T[m] >= t).sum())}", ha="center", va="center", fontsize=7.2,
                        transform=ax.get_xaxis_transform(), color=col)
        lr = logrank_test(T[hi.values == 1], T[hi.values == 0], E[hi.values == 1], E[hi.values == 0])
        cph = CoxPHFitter().fit(pd.DataFrame({"T": T, "E": E, "g": hi.values}), "T", "E")
        hr = np.exp(cph.params_["g"]); lo_, up_ = np.exp(cph.confidence_intervals_.values[0])
        res.append(dict(cohort=c, method=meth, n_high=int(hi.sum()), hr=hr, lo=lo_, up=up_, p=lr.p_value))
        ptxt = "P<0.001" if lr.p_value < 0.001 else (f"P={lr.p_value:.3f}" if lr.p_value < 0.01 else f"P={lr.p_value:.2f}")
        ax.text(0.03, 0.05, f"HR {hr:.2f} (95% CI {lo_:.2f}–{up_:.2f})\nlog-rank {ptxt}", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=7)
        ax.set_xlim(0, 62); ax.set_xticks(ticks); ax.set_ylim(0, 1.02)
        ax.set_xlabel("Months", labelpad=1); ax.set_title(f"{'Gastric' if c == 'STAD' else 'Colorectal'} cancer, {lab}", fontsize=8)
        ax.get_legend().remove()
        ax.text(-0.05, -0.36, "Low", ha="right", va="center", fontsize=7.2, transform=ax.transAxes, color=LO)
        ax.text(-0.05, -0.45, "High", ha="right", va="center", fontsize=7.2, transform=ax.transAxes, color=HI)
        ax.text(-0.02, -0.27, "No. at risk", ha="right", va="center", fontsize=7.2, transform=ax.transAxes)
        ax.text(-0.13, 1.06, "ABCD"[2 * r + k], transform=ax.transAxes, fontsize=8, fontweight="bold")
    axes[r, 0].set_ylabel("Overall survival")
from matplotlib.lines import Line2D
fig.legend([Line2D([], [], color=LO, ls="--"), Line2D([], [], color=HI)], ["Low predicted risk", "High predicted risk"],
           loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.savefig(str(OUT / "study2_km_median_split.png"), dpi=600); fig.savefig(str(OUT / "study2_km_median_split.tif"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
R = pd.DataFrame(res); print(R.round(3).to_string())
from sklearn.metrics import cohen_kappa_score
for c in ["STAD", "CRC"]:
    a, b = groups[(c, "Clinical-Cox")], groups[(c, "Cox+DCA")]
    print(c, "agreement %.1f%%" % (100 * (a == b).mean()), "kappa %.2f" % cohen_kappa_score(a, b), "n", len(a))
R.to_csv(REPO / "results" / "study2" / "km_median_split.csv", index=False)
