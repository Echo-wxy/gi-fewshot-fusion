# Experiment 2, protocol v2 — STAD (encoder uni, tiles <= 0)
Patients 278, deaths 108, censoring 61%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.646 ± 0.049 | 0.640 | 0.174 | 0.162 | 0.220 | 0.649 |
| Cox+zero-residual | 0.646 ± 0.049 | 0.640 | 0.187 | 0.178 | 0.229 | 0.649 |
| Omics-CoxNet | 0.547 ± 0.062 | 0.531 | 0.195 | 0.171 | 0.258 | 0.552 |
| Stacked-Cox | 0.638 ± 0.033 | 0.634 | 0.175 | 0.163 | 0.222 | 0.646 |
| Late-Fusion | 0.612 ± 0.047 | 0.601 | — | — | — | 0.627 |
| WSI-only | 0.530 ± 0.060 | 0.525 | 0.200 | 0.183 | 0.253 | 0.540 |
| Omics-only | 0.546 ± 0.062 | 0.535 | 0.202 | 0.183 | 0.258 | 0.548 |
| Concat-Fusion | 0.578 ± 0.065 | 0.574 | 0.194 | 0.179 | 0.243 | 0.587 |
| Attention-Fusion | 0.550 ± 0.072 | 0.542 | 0.195 | 0.178 | 0.244 | 0.568 |
| Gated-Fusion | 0.564 ± 0.069 | 0.565 | 0.196 | 0.179 | 0.246 | 0.582 |
| DCA | 0.559 ± 0.076 | 0.556 | 0.195 | 0.178 | 0.246 | 0.575 |
| PORPOISE-style | 0.568 ± 0.064 | 0.566 | 0.191 | 0.173 | 0.250 | 0.565 |
| MCAT-style | 0.543 ± 0.070 | 0.527 | 0.194 | 0.176 | 0.250 | 0.551 |
| Cox+WSI | 0.630 ± 0.042 | 0.625 | 0.189 | 0.180 | 0.229 | 0.638 |
| Cox+Omics | 0.647 ± 0.043 | 0.636 | 0.188 | 0.179 | 0.230 | 0.654 |
| Cox+Concat | 0.635 ± 0.049 | 0.630 | 0.190 | 0.180 | 0.232 | 0.641 |
| Cox+Attention | 0.626 ± 0.043 | 0.618 | 0.189 | 0.178 | 0.230 | 0.639 |
| Cox+Gated | 0.631 ± 0.054 | 0.626 | 0.191 | 0.179 | 0.233 | 0.647 |
| Cox+DCA | 0.633 ± 0.050 | 0.629 | 0.189 | 0.179 | 0.230 | 0.648 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | -0.001 | [-0.034, +0.034] | 0.517 |
| Cox+DCA − Cox+Attention | +0.009 | [-0.008, +0.026] | 0.136 |
| Cox+DCA − Cox+Concat | +0.007 | [-0.010, +0.025] | 0.194 |
| Cox+DCA − Cox+Gated | +0.001 | [-0.010, +0.013] | 0.461 |
| Cox+DCA − Stacked-Cox | +0.002 | [-0.028, +0.033] | 0.453 |
| Cox+DCA − Late-Fusion | +0.021 | [-0.009, +0.056] | 0.086 |
| Cox+WSI − Clinical-Cox | -0.011 | [-0.040, +0.017] | 0.769 |
| Cox+Omics − Clinical-Cox | +0.005 | [-0.029, +0.042] | 0.395 |
| DCA − Attention-Fusion | +0.007 | [-0.017, +0.032] | 0.280 |
| DCA − Concat-Fusion | -0.011 | [-0.040, +0.019] | 0.770 |
| DCA − Gated-Fusion | -0.007 | [-0.025, +0.010] | 0.783 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 12 | 32 | 1.117 | 0.826 | 1.06 | 0.6 |
| Concat-Fusion | 3 | 23 | 1.135 | 0.658 | 1.78 | 0.4 |
| Cox+Attention | 7 | 27 | 1.096 | 0.886 | 1.09 | 0.5 |
| Cox+Concat | 2 | 22 | 1.119 | 0.683 | 1.77 | 0.4 |
| Cox+DCA | 9 | 29 | 1.117 | 0.758 | 1.08 | 0.8 |
| Cox+Gated | 6 | 26 | 1.115 | 0.655 | 1.15 | 0.6 |
| Cox+Omics | 2 | 22 | 1.090 | 0.996 | 0.69 | 0.3 |
| Cox+WSI | 2 | 22 | 1.113 | 1.008 | 0.59 | 0.3 |
| Cox+zero-residual | 98 | 100 | 1.091 | 1.106 | 0.17 | 1.3 |
| DCA | 13 | 33 | 1.126 | 0.729 | 1.04 | 0.8 |
| Gated-Fusion | 8 | 28 | 1.143 | 0.641 | 1.15 | 0.6 |
| MCAT-style | 7 | 27 | 1.131 | 0.543 | 2.26 | 6.7 |
| PORPOISE-style | 8 | 28 | 1.147 | 0.548 | 2.46 | 6.6 |
| stage1-clin | 98 | 100 | 1.138 | 1.110 | 0.39 | 1.7 |
| stage1-omics | 13 | 33 | 1.154 | 1.015 | 2.00 | 0.5 |
| stage1-wsi | 2 | 22 | 1.169 | 0.713 | 2.11 | 4.9 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.364 | 0.313 | 0.322 | 0.871 | -0.034 | -0.005 | 0.001 | 3.708 | 3.949 | 1.705 |
| Cox+Attention | 0.355 | 0.320 | 0.325 | 0.871 | -0.035 | 0.002 | -0.004 | 3.708 | 3.949 | 1.705 |
| Cox+DCA | 0.377 | 0.305 | 0.318 | 0.887 | -0.054 | -0.016 | -0.005 | 3.708 | 3.949 | 1.705 |
| DCA | 0.374 | 0.307 | 0.319 | 0.883 | -0.055 | -0.011 | 0.005 | 3.708 | 3.949 | 1.705 |