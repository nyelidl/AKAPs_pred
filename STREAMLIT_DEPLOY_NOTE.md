# Streamlit Cloud deployment note

1. Upload/commit all files in this folder to GitHub.
2. On Streamlit Cloud, select this repository and set the main file to:

```text
akap_app.py
```

3. The app expects these files in the repository root:

```text
akap_screen.py
akap_ml.py
akap_pssm.json
akap_ml_model_v2.joblib
akap_ml_model_v2_metadata.json
```

A compatibility copy named `akap_ml_model.joblib` is also included.

4. This is AKAPSpred v5.1 using ML v2. ML v3 is not deployed yet.
