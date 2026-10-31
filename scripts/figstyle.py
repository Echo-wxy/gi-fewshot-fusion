import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "Liberation Serif", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 1.2,
    "savefig.dpi": 600, "pdf.fonttype": 42})
MM = 1 / 25.4
# color palette (one color, one meaning)
TEAL, TEAL_E = "#2E9E86", "#14826B"     # proposed adapter (Study 1)
BLUE, BLUE_E = "#4C7DC0", "#2F62A8"     # clinical information / clinical Cox
PURPLE = "#8E6FB8"                       # WSI / pathology
AMBER = "#E3A03A"                        # omics
GRAY = "#8A8A8A"                         # compared / training-free baselines
RED = "#C8504A"                          # negative finding
