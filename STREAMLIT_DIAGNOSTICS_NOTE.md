# Streamlit diagnostics note — confirmed danger-zone false positive

_For the Model Diagnostics view (or release notes). AKAPSpred v5.1 + ML v2, frozen._

**Known false positive (real protein, logged 2026-06-12).** A screen of
phosphodiesterases produced one high-confidence call that is a contextual false
positive:

- **PDE2A 404–427** (core `VSVLLQEIITEA`), isoform RII, PSSM 15.14, ML 0.992 →
  promoted to `high` via the high-background downgrade (pssm ≥ 14 AND ml ≥ 0.95).
- PDE2A is a phosphodiesterase; this window lies in the **GAF-B regulatory /
  dimerization domain**, not a PKA RI/RII D/D anchoring helix. Classified as a
  *biologically supported non-AKAP contextual false positive* (domain/function
  annotation) — **not** an experimentally proven non-binder; no direct PKA binding
  assay exists, and the burden of proof is on the AKAP-positive call.

**Correctly rejected (same screen):** PDE3A 64–87 (sensitive_only; rejected by all
three layers), PDE4D 267–290 (sensitive_only), PDE2A 302–325 (unlikely,
negative-determinant flag). So this does **not** mean all PDE helices fail.

**What it means for users.** A `high` call in the **danger zone**
(`amphipathic = True` and `PSSM ≥ 12`) is not yet independently checked by ML —
ML v2 can endorse a real amphipathic non-AKAP there (here, 0.99). Treat
high-confidence hits in that region, especially in regulatory/dimerization or
GAF-domain contexts, as hypotheses for experimental validation. This is the
priority target for the planned ML v3 hard negatives; v5.1 + ML v2 remains the
deployed model.
