#!/bin/sh
# Runs the validation and all figure scripts. Figures are saved in the current folder.
set -e
export MPLBACKEND=Agg
python3 validate_wkb.py
python3 figure1_electron_proton.py
python3 figure2_enzyme.py
python3 figure3_KIE.py
python3 figure4_arrhenius_final.py
python3 figure5.py      # about 20 minutes
python3 figure6_zpe.py
mkdir -p figures && mv -f *.png figures/ 2>/dev/null || true
