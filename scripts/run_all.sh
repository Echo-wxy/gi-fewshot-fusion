#!/usr/bin/env bash
# Full pipeline. Requires a CUDA GPU, about 60 GB of disk in $GIFS_ROOT and HF_TOKEN with access to MahmoodLab/UNI.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${GIFS_ROOT:?set GIFS_ROOT to a working directory}"

python src/download_data.py

# Study 1
python src/exp1_endo.py all
python src/exp1_endo.py run --mask-overlay
python src/exp1_endo.py figures --mask-overlay

# Study 2
python src/exp2_extract.py --encoders uni,resnet50,phikon
python src/exp2_model.py prep
python src/exp2_model.py taskA
python src/exp2_v2.py --encoder uni --cap 500
python src/exp2_v2.py --encoder resnet50 --cap 500
python src/exp2_v2.py --encoder phikon --cap 500
python src/exp2_v2.py --encoder uni --cap 100
python src/exp2_v2.py --encoder uni --cap 0
