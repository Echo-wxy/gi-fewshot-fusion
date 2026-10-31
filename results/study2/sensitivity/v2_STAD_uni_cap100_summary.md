# Experiment 2, protocol v2 — STAD (encoder uni, tiles <= 100)
Patients 278, deaths 108, censoring 61%; outer 5-fold x 3 seeds; inner 3-fold ensembles.

| Method | Harrell C (mean ± SD) | Uno C | IBS | Brier 12 m | Brier 24 m | Pooled C |
|---|---|---|---|---|---|---|
| Clinical-Cox | 0.646 ± 0.049 | 0.640 | 0.174 | 0.162 | 0.220 | 0.649 |
| Cox+zero-residual | 0.646 ± 0.049 | 0.640 | 0.187 | 0.178 | 0.229 | 0.649 |
| Omics-CoxNet | 0.547 ± 0.062 | 0.531 | 0.195 | 0.171 | 0.258 | 0.552 |
| Stacked-Cox | 0.639 ± 0.031 | 0.634 | 0.175 | 0.163 | 0.222 | 0.646 |
| Late-Fusion | 0.615 ± 0.047 | 0.603 | — | — | — | 0.630 |
| WSI-only | 0.533 ± 0.061 | 0.530 | 0.199 | 0.184 | 0.252 | 0.544 |
| Omics-only | 0.546 ± 0.062 | 0.535 | 0.201 | 0.183 | 0.258 | 0.549 |
| Concat-Fusion | 0.577 ± 0.068 | 0.573 | 0.195 | 0.179 | 0.244 | 0.587 |
| Attention-Fusion | 0.547 ± 0.074 | 0.540 | 0.196 | 0.178 | 0.245 | 0.564 |
| Gated-Fusion | 0.561 ± 0.070 | 0.560 | 0.197 | 0.180 | 0.248 | 0.577 |
| DCA | 0.559 ± 0.072 | 0.554 | 0.195 | 0.178 | 0.247 | 0.574 |
| PORPOISE-style | 0.572 ± 0.067 | 0.571 | 0.191 | 0.172 | 0.250 | 0.572 |
| MCAT-style | 0.542 ± 0.069 | 0.527 | 0.194 | 0.176 | 0.250 | 0.551 |
| Cox+WSI | 0.632 ± 0.041 | 0.627 | 0.189 | 0.180 | 0.228 | 0.641 |
| Cox+Omics | 0.648 ± 0.043 | 0.637 | 0.188 | 0.178 | 0.230 | 0.655 |
| Cox+Concat | 0.630 ± 0.041 | 0.624 | 0.191 | 0.179 | 0.233 | 0.637 |
| Cox+Attention | 0.624 ± 0.044 | 0.616 | 0.189 | 0.178 | 0.230 | 0.636 |
| Cox+Gated | 0.633 ± 0.053 | 0.627 | 0.191 | 0.179 | 0.233 | 0.649 |
| Cox+DCA | 0.633 ± 0.048 | 0.629 | 0.189 | 0.179 | 0.231 | 0.649 |

Primary contrast (pre-specified): Cox+DCA − Clinical-Cox, pooled out-of-fold C, patient-level bootstrap.

| Contrast | ΔC | 95% CI | Share of bootstrap ≤ 0 |
|---|---|---|---|
| Cox+DCA − Clinical-Cox | -0.000 | [-0.033, +0.033] | 0.500 |
| Cox+DCA − Cox+Attention | +0.012 | [-0.004, +0.028] | 0.067 |
| Cox+DCA − Cox+Concat | +0.011 | [-0.004, +0.028] | 0.085 |
| Cox+DCA − Cox+Gated | -0.001 | [-0.011, +0.010] | 0.570 |
| Cox+DCA − Stacked-Cox | +0.003 | [-0.026, +0.033] | 0.432 |
| Cox+DCA − Late-Fusion | +0.018 | [-0.012, +0.053] | 0.121 |
| Cox+WSI − Clinical-Cox | -0.008 | [-0.037, +0.019] | 0.704 |
| Cox+Omics − Clinical-Cox | +0.006 | [-0.029, +0.043] | 0.382 |
| DCA − Attention-Fusion | +0.010 | [-0.015, +0.036] | 0.217 |
| DCA − Concat-Fusion | -0.013 | [-0.043, +0.019] | 0.778 |
| DCA − Gated-Fusion | -0.003 | [-0.019, +0.013] | 0.639 |

Training diagnostics (median over all inner models):

| Model | Best epoch | Stop epoch | Best val NLL | Last train NLL | Grad norm | Seconds |
|---|---|---|---|---|---|---|
| Attention-Fusion | 11 | 31 | 1.119 | 0.834 | 1.06 | 0.5 |
| Concat-Fusion | 3 | 23 | 1.141 | 0.663 | 1.77 | 0.4 |
| Cox+Attention | 7 | 27 | 1.094 | 0.884 | 1.08 | 0.5 |
| Cox+Concat | 2 | 22 | 1.115 | 0.683 | 1.78 | 0.4 |
| Cox+DCA | 10 | 30 | 1.119 | 0.752 | 1.09 | 0.7 |
| Cox+Gated | 6 | 26 | 1.113 | 0.658 | 1.15 | 0.6 |
| Cox+Omics | 2 | 22 | 1.090 | 0.991 | 0.69 | 0.3 |
| Cox+WSI | 3 | 23 | 1.112 | 1.009 | 0.59 | 0.3 |
| Cox+zero-residual | 98 | 100 | 1.091 | 1.106 | 0.17 | 1.2 |
| DCA | 12 | 32 | 1.127 | 0.721 | 1.05 | 0.8 |
| Gated-Fusion | 8 | 28 | 1.138 | 0.638 | 1.17 | 0.6 |
| MCAT-style | 7 | 27 | 1.138 | 0.496 | 2.25 | 1.4 |
| PORPOISE-style | 7 | 27 | 1.143 | 0.560 | 2.37 | 1.2 |
| stage1-clin | 98 | 100 | 1.138 | 1.110 | 0.39 | 1.6 |
| stage1-omics | 13 | 33 | 1.154 | 1.013 | 2.00 | 0.5 |
| stage1-wsi | 2 | 22 | 1.167 | 0.710 | 2.13 | 0.7 |

Fusion mechanism (mean over folds; dNLL_remove_* = rise in test NLL when that modality is replaced by its training mean):

| Method | alpha_wsi | alpha_omics | alpha_clin | entropy | dNLL_remove_wsi | dNLL_remove_omics | dNLL_remove_clin | norm_wsi | norm_omics | norm_clin |
|---|---|---|---|---|---|---|---|---|---|---|
| Attention-Fusion | 0.374 | 0.308 | 0.318 | 0.870 | -0.031 | 0.002 | 0.001 | 3.712 | 3.945 | 1.704 |
| Cox+Attention | 0.360 | 0.317 | 0.322 | 0.875 | -0.029 | 0.005 | -0.003 | 3.712 | 3.945 | 1.704 |
| Cox+DCA | 0.377 | 0.309 | 0.314 | 0.889 | -0.037 | 0.007 | -0.001 | 3.712 | 3.945 | 1.704 |
| DCA | 0.375 | 0.305 | 0.321 | 0.878 | -0.062 | -0.015 | 0.004 | 3.712 | 3.945 | 1.704 |