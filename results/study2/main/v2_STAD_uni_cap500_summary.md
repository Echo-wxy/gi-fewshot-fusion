# Experiment 2, protocol v2 — STAD (encoder uni, tiles <= 500)
Patients 278, deaths 108, censoring 61%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.646 ± 0.049 | 0.640 | 0.174 | 0.162 | 0.220 | 0.649 |
| Cox+zero-residual | 0.646 ± 0.049 | 0.640 | 0.187 | 0.178 | 0.229 | 0.649 |
| Omics-CoxNet | 0.547 ± 0.062 | 0.531 | 0.195 | 0.171 | 0.258 | 0.552 |
| Stacked-Cox | 0.638 ± 0.031 | 0.634 | 0.175 | 0.163 | 0.222 | 0.646 |
| Late-Fusion | 0.612 ± 0.048 | 0.601 | — | — | — | 0.627 |
| WSI-only | 0.530 ± 0.061 | 0.526 | 0.200 | 0.184 | 0.253 | 0.541 |
| Omics-only | 0.547 ± 0.062 | 0.535 | 0.201 | 0.183 | 0.258 | 0.550 |
| Concat-Fusion | 0.581 ± 0.062 | 0.578 | 0.194 | 0.179 | 0.243 | 0.592 |
| Attention-Fusion | 0.550 ± 0.072 | 0.541 | 0.195 | 0.178 | 0.245 | 0.569 |
| Gated-Fusion | 0.564 ± 0.066 | 0.564 | 0.196 | 0.179 | 0.247 | 0.583 |
| DCA | 0.558 ± 0.073 | 0.556 | 0.196 | 0.179 | 0.248 | 0.575 |
| PORPOISE-style | 0.570 ± 0.066 | 0.568 | 0.192 | 0.173 | 0.250 | 0.568 |
| MCAT-style | 0.547 ± 0.071 | 0.532 | 0.193 | 0.175 | 0.248 | 0.556 |
| Cox+WSI | 0.630 ± 0.042 | 0.625 | 0.189 | 0.180 | 0.229 | 0.638 |
| Cox+Omics | 0.648 ± 0.043 | 0.637 | 0.188 | 0.178 | 0.230 | 0.656 |
| Cox+Concat | 0.631 ± 0.042 | 0.626 | 0.190 | 0.179 | 0.232 | 0.639 |
| Cox+Attention | 0.626 ± 0.043 | 0.618 | 0.188 | 0.178 | 0.229 | 0.639 |
| Cox+Gated | 0.631 ± 0.052 | 0.625 | 0.191 | 0.179 | 0.233 | 0.646 |
| Cox+DCA | 0.630 ± 0.052 | 0.627 | 0.190 | 0.179 | 0.232 | 0.646 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | -0.003 | [-0.036, +0.033] | 0.551 |
| Cox+DCA − Cox+Attention | +0.008 | [-0.009, +0.024] | 0.178 |
| Cox+DCA − Cox+Concat | +0.007 | [-0.009, +0.024] | 0.174 |
| Cox+DCA − Cox+Gated | -0.000 | [-0.011, +0.012] | 0.518 |
| Cox+DCA − Stacked-Cox | +0.001 | [-0.030, +0.032] | 0.474 |
| Cox+DCA − Late-Fusion | +0.019 | [-0.010, +0.053] | 0.105 |
| Cox+WSI − Clinical-Cox | -0.010 | [-0.040, +0.018] | 0.751 |
| Cox+Omics − Clinical-Cox | +0.007 | [-0.028, +0.043] | 0.368 |
| DCA − Attention-Fusion | +0.006 | [-0.018, +0.030] | 0.324 |
| DCA − Concat-Fusion | -0.017 | [-0.048, +0.015] | 0.862 |
| DCA − Gated-Fusion | -0.008 | [-0.027, +0.009] | 0.810 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 12 | 32 | 1.123 | 0.834 | 1.06 | 0.6 |
| Concat-Fusion | 3 | 23 | 1.134 | 0.665 | 1.77 | 0.4 |
| Cox+Attention | 7 | 27 | 1.096 | 0.888 | 1.08 | 0.5 |
| Cox+Concat | 3 | 23 | 1.116 | 0.689 | 1.77 | 0.4 |
| Cox+DCA | 11 | 31 | 1.123 | 0.762 | 1.07 | 0.8 |
| Cox+Gated | 6 | 26 | 1.119 | 0.658 | 1.15 | 0.5 |
| Cox+Omics | 2 | 22 | 1.090 | 0.991 | 0.69 | 0.3 |
| Cox+WSI | 3 | 23 | 1.111 | 1.001 | 0.58 | 0.3 |
| Cox+zero-residual | 98 | 100 | 1.091 | 1.106 | 0.17 | 1.2 |
| DCA | 13 | 33 | 1.138 | 0.733 | 1.04 | 0.8 |
| Gated-Fusion | 8 | 28 | 1.133 | 0.644 | 1.13 | 0.6 |
| MCAT-style | 7 | 27 | 1.136 | 0.516 | 2.25 | 2.3 |
| PORPOISE-style | 7 | 27 | 1.139 | 0.539 | 2.46 | 2.1 |
| stage1-clin | 98 | 100 | 1.138 | 1.110 | 0.39 | 1.6 |
| stage1-omics | 13 | 33 | 1.154 | 1.015 | 2.00 | 0.6 |
| stage1-wsi | 2 | 22 | 1.168 | 0.712 | 2.12 | 1.5 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.357 | 0.315 | 0.327 | 0.873 | -0.031 | -0.002 | 0.001 | 3.690 | 3.949 | 1.705 |
| Cox+Attention | 0.350 | 0.320 | 0.330 | 0.876 | -0.034 | 0.005 | -0.003 | 3.690 | 3.949 | 1.705 |
| Cox+DCA | 0.379 | 0.309 | 0.312 | 0.886 | -0.051 | -0.006 | -0.002 | 3.690 | 3.949 | 1.705 |
| DCA | 0.371 | 0.308 | 0.321 | 0.886 | -0.053 | -0.012 | 0.003 | 3.690 | 3.949 | 1.705 |