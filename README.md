# gi-fewshot-fusion

Code and result files for two analyses on public gastrointestinal data:

1. **Few-shot recognition of rare lesions in endoscopic images (HyperKvasir).** A zero-initialized residual adapter on frozen DINOv2 features is trained on common classes and evaluated on held-out rare lesion classes in 1-, 5- and 10-shot episodes, with near-duplicate frames kept apart.
2. **Survival prediction from whole-slide image (WSI) tumor tiles, omics and clinical data (TCGA gastric and colorectal cancer).** Nineteen survival models, including differentiable conflict-aware attention (DCA) fusion anchored on a clinical Cox model, are compared with the clinical Cox model under nested, repeated cross-validation.

## Repository layout

```
src/
  download_data.py   resumable download of all public data
  exp1_endo.py       Study 1: feature extraction, near-duplicate grouping, episodic evaluation, figures
  exp2_extract.py    Study 2: tile features with UNI, ResNet50 or Phikon
  exp2_model.py      Study 2: data preparation, encoders, metrics, MSI/hypermutation tile check
  exp2_v2.py         Study 2: final survival protocol (clinical offset, two-stage fusion, diagnostics)
scripts/
  run_all.sh                 full pipeline
  make_tables.py             main tables from results/
  per_class.py               per-class 5-shot accuracy from the query-level predictions
  plot_*.py                  figures from results/
results/
  study1/   summaries, paired differences, splits, near-duplicate clusters, query-level predictions
  study2/   per-setting summaries, contrasts, fold metrics, diagnostics, fusion weights,
            out-of-fold predictions, cohort tables, MSI tile check, median-split survival analysis
tests/      consistency checks on the result files
data/       data sources and licenses
```

## Quick start (no GPU, no raw data)

```bash
pip install -r requirements.txt
python scripts/make_tables.py
python scripts/per_class.py
python scripts/plot_study1.py
python scripts/plot_contrasts.py
python scripts/plot_fusion_mechanism.py
python scripts/plot_km.py
python -m pytest -q tests
```

Figures are written to `figures/`. All commands finish in under a minute on a laptop.

## Full pipeline

Hardware used: one NVIDIA RTX 4070 (12 GB). Total run time is about one day, most of it tile feature extraction and the Study 2 sensitivity runs.

```bash
export GIFS_ROOT=/path/to/workdir        # data, features and results are written here (about 60 GB)
export HF_TOKEN=...                      # Hugging Face token with access to MahmoodLab/UNI
bash scripts/run_all.sh
```

`run_all.sh` runs, in order: data download, Study 1 (`exp1_endo.py all`, then the overlay-masked repeat), tile feature extraction, Study 2 data preparation and MSI tile check (`exp2_model.py prep` and `taskA`), the main survival analysis (`exp2_v2.py`) and the four sensitivity runs (ResNet50, Phikon, at most 100 tiles, all tiles).

UNI is a gated model: request access on its Hugging Face page before running the extraction step.

## Main results (values computed by `scripts/make_tables.py` from `results/`)

| Analysis | Value |
|---|---|
| Study 1, 5-shot accuracy on five rare lesion classes, episodic adapter | 61.5% |
| Study 1, same setting, SimpleShot | 55.3% |
| Study 2, pooled C-index of the clinical Cox model, gastric / colorectal | 0.649 / 0.687 |
| Study 2, difference Cox+DCA minus clinical Cox, gastric | −0.003 (95% CI −0.036 to 0.033) |
| Study 2, difference Cox+DCA minus clinical Cox, colorectal | 0.018 (95% CI −0.032 to 0.072) |

## Method names in the result files

| File label | Description |
|---|---|
| `ProtoNet-Head (ours)` | episodically trained residual adapter |
| `NonEpisodic-Head+NC` | same adapter trained with class-balanced cross-entropy |
| `SimpleShot-CL2N`, `NC-L2N`, `LogReg-CL2N` | training-free baselines on the same features |
| `Clinical-Cox` | Cox model on age, sex and pathological stage |
| `Cox+...` | network with the clinical Cox linear predictor as offset |
| `DCA` | differentiable conflict-aware attention fusion |
| `PORPOISE-style`, `MCAT-style` | re-implementations of published architectures in this framework |

## Data

All inputs are public. Sources and licenses are listed in `data/README.md`.

## License

Code: MIT (see `LICENSE`). Result files: CC BY 4.0. Third-party data keep their original licenses.
