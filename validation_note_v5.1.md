# AKAPSpred v5.1 validation note

AKAPSpred v5.1 separates sensitive PSSM-based motif discovery from confidence-ranked proteomic prioritization. The v5.1 confidence logic restores high-confidence recovery of synthetic RI/RII AKAP-like positives while keeping high-confidence calls in random, composition-shuffled, and DDIP-like backgrounds near zero.

Selected thresholds:

- RI_PSSM_HIGH_THR = 12.0
- RII_PSSM_HIGH_THR = 12.0
- ML_HIGH_THR = 0.80
- ML_VHIGH_THR = 0.90
- ML_MEDIUM_THR = 0.60
- ML_HIGH_BG_THR = 0.95
- HIGH_BG_PSSM_BONUS = 2.0

High/elevated background risk is treated as a downgrade rather than a hard block: such hits can reach `high` only with stricter ML support and stronger PSSM support, and they cannot reach `very_high`.

Important limitation: the negative_determinant_disrupted synthetic class is not fully resolved. Window-level validation shows that many residual high-confidence cases are on-target and arise from PSSM register/frame offset, where the introduced Asp/Glu falls outside the exact anchor frame checked by the determinant detector. Broad register-tolerant determinant scanning was rejected because it falsely flags many true positives. This limitation is reported explicitly rather than hidden with an overfitted rule.

This synthetic benchmark is for internal development and stress-testing; it is not external biological validation.
