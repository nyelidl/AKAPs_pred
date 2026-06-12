# AKAPSpred v5.1 — confidence-logic changes

## Summary
Fixes over-strict final-confidence logic that suppressed true RII-like positives,
adds isoform-specific thresholds, adds window-level on/off-target validation,
and reruns the synthetic 10k benchmark. PSSM sensitivity unchanged; ML model unchanged.

## akap_screen.py
- `background_risk` of "elevated"/"high" is no longer a HARD BLOCK. It is now a
  DOWNGRADE path: such hits can reach `high` (never `very_high`) but only via a
  stricter ML gate (ML_HIGH_BG_THR) AND a raised PSSM floor (HIGH_BG_PSSM_BONUS).
- Low/moderate background keeps original semantics: very_high (ml>=0.90) / high (ml>=0.80).
- Isoform-specific PSSM thresholds via `_iso_high_thr(iso)`.
- Red flags (unlikely class, n_negdet>=2, or n_negdet==1 with weak PSSM) still force `unlikely`;
  ML cannot override them.

### Shipped thresholds
RII_PSSM_HIGH_THR=12.0, RI_PSSM_HIGH_THR=12.0,
ML_HIGH_THR=0.80, ML_VHIGH_THR=0.90, ML_MEDIUM_THR=0.60,
ML_HIGH_BG_THR=0.95, HIGH_BG_PSSM_BONUS=2.0

## validate_synthetic_benchmark_v2.py
- Window-level overlap vs embedded motif (motif_start/motif_end): overlap_length,
  overlap_fraction, overlaps_embedded_motif, on_target_hit, off_target_hit.
- Per-category on/off-target any-hit and high-confidence rates.
- New outputs: window_level_summary.csv, predictions_with_window_overlap.csv.
- Backward-compatible (old outputs still produced).

## tune_confidence_thresholds.py (new)
- Re-derives tiers in-memory from a baseline prediction CSV across a 288-cell grid.
- Outputs threshold_tuning_summary.csv and the recommended set.

## Notebook
- AKAPSpred_synthetic_10k_validation_colab_v51.ipynb: v5.1 note + window_level_summary display.

## Streamlit diagnostics tab (this update)
- akap_app.py: added a sidebar "View" toggle and a read-only "📊 Model diagnostics"
  page with: deployed-model card (v5.1 / ML v2 / v3 not deployed), the PDE sanity-check
  table, the PDE2A 404–427 contextual false-positive note (annotation-based, not an
  experimental non-binder; burden of proof on the AKAP-positive call), the danger-zone
  limitation (amphipathic & PSSM≥12 & non-AKAP context), the ML v3 hard-negative
  direction, and user guidance (high-confidence = candidate, review danger-zone hits).
- Prediction behavior, confidence logic, thresholds, and the ML v2 model are UNCHANGED.
  The diagnostics view renders then st.stop()s; the screener path is not touched.
- No PDE-specific rejection rule was added.
