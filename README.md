# AKAPSpred v5.1

AKAPSpred v5.1 screens protein sequences for AKAP-like PKA regulatory-subunit anchoring helices.

This repository package uses **AKAPSpred v5.1 confidence logic with ML v2**.

## Main files

- `akap_app.py` — Streamlit web app
- `akap_screen.py` — command-line/backend screener
- `akap_from_csv.py` — CSV/batch helper
- `akap_ml.py` — ML feature and model utilities
- `akap_pssm.json` — RI/RII PSSM profiles
- `akap_ml_model_v2.joblib` — current ML v2 model
- `akap_ml_model.joblib` — compatibility copy of the v2 model
- `validate_akap_screen.py` — Burgers/THAHIT SI validation helper
- `validate_synthetic_benchmark_v2.py` — synthetic benchmark validation with window-level analysis

## Run Streamlit

```bash
pip install -r requirements.txt
streamlit run akap_app.py
```

## Run command-line screening

```bash
python akap_screen.py input.fasta --use-ml --ml-model akap_ml_model_v2.joblib -o predictions.csv
```

## Important note

The deployed ML model is **ML v2**. AKAPSpred v5.1 improves confidence logic, validation, and documentation but does not replace the ML engine. ML probabilities should be interpreted cautiously as ML confidence/ranking scores unless calibration diagnostics support probability interpretation.

See:

- `CHANGES_v5.1.md`
- `validation_note_v5.1.md`
