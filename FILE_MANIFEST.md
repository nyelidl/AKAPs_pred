# AKAPs_pred v5.1 GitHub/Streamlit package

This package is built from the current v5.1 repo files plus the newly uploaded biological-context hard-negative benchmark files from `files_new.zip`.

## Core app/repo files
- README.md
- CHANGES_v5.1.md
- STREAMLIT_DEPLOY_NOTE.md
- akap_app.py
- akap_from_csv.py
- akap_ml.py
- akap_ml_model.joblib
- akap_ml_model_v2.joblib
- akap_ml_model_v2_metadata.json
- akap_pssm.json
- akap_screen.py
- index.html
- requirements.txt
- train_akap_ml_v2.py
- tune_confidence_thresholds.py
- validate_akap_screen.py
- validate_synthetic_benchmark_v2.py
- validation_note_v5.1.md
- before_after_10k_validation.csv
- threshold_tuning_summary.csv
- window_level_summary.csv

## Newly uploaded biological-context benchmark files
- BIOLOGICAL_HARD_NEGATIVE_BENCHMARK_PROTOCOL.md
- hard_negative_candidate_list_biological_context.md
- hard_negative_candidate_list_biological_context.csv
- hard_negative_metadata_template_biological_context.csv
- run_hard_negative_benchmark.py
- hard_negative_qc.py
- hard_negative_readme.md
- build_candidate_list.py

## Version note
Current deployed model remains AKAPSpred v5.1 + ML v2. The biological hard-negative benchmark files are diagnostic/ML-v3 preparation files and do not change prediction behavior.
