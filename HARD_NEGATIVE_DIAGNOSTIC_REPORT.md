# Hard-negative diagnostic report — AKAPSpred v5.1 + ML v2

_Logged from a real screen of phosphodiesterases (user-supplied export, 2026-06-12). Model frozen; no retraining, no threshold change._

Danger zone = `amphipathic == True AND pssm_score >= 12`. PDEs are cyclic-nucleotide effectors in cAMP/PKA compartments but are **not** AKAPs, so any high/very_high call here is a contextual false positive (annotation-based).

## Windows screened

| protein      |   start |   end | core         | isoform   |   pssm_score |   ml_prob | amphipathic   |   n_negdet | proteomic_confidence   | danger   | is_fp   |
|:-------------|--------:|------:|:-------------|:----------|-------------:|----------:|:--------------|-----------:|:-----------------------|:---------|:--------|
| PDE2A_O00408 |     404 |   427 | VSVLLQEIITEA | RII       |        15.14 |    0.9917 | True          |          0 | high                   | True     | True    |
| PDE4D_Q08499 |     267 |   290 | YQKLASETLEEL | RII       |        11.04 |    0.2715 | True          |          0 | sensitive_only         | False    | False   |
| PDE2A_O00408 |     302 |   325 | LKDLTSEDVQQL | RII       |         8.15 |    0.0851 | True          |          1 | unlikely               | False    | False   |
| PDE3A_Q14432 |      64 |    87 | LSFLLALLVRLV | RII       |         7.64 |    0.005  | False         |          0 | sensitive_only         | False    | False   |
| PDE2A_O00408 |      57 |    80 | LQRAVKEALSAV | RII       |         7.4  |    0.0078 | True          |          0 | sensitive_only         | False    | False   |

## Class: PDE_GAF_domain_dimerization_helix

- proteins: 3 (PDE2A, PDE3A, PDE4D)
- windows screened: 5
- danger-zone windows: 1
- danger-zone ML median: 0.992; fraction ML>=0.95: 1.00
- high-confidence false positives: 1 (all 1 in the danger zone, full-protein, in known helix region)
- **verdict: PRIORITY v3 hard-negative class — enters the danger zone AND ML v2 endorses it (0.99).**

## Key interpretation

- **PDE3A 64-87 is NOT a false positive.** It is rejected by all three layers (PSSM 7.64 < 12, ML 0.005, amphipathic = False) -> sensitive_only.
- **PDE2A 404-427 (VSVLLQEIITEA) is the important finding:** a real full-protein, high-PSSM (15.14), amphipathic non-AKAP that ML v2 scored 0.992 and the pipeline promoted to `high` via the high-background downgrade (pssm>=14 AND ml>=0.95). This is a **biologically supported non-AKAP contextual false positive** (the window lies in the GAF-B regulatory/dimerization domain, not a PKA RI/RII D/D anchoring helix). It is NOT an experimentally proven non-binder — no direct PKA binding assay exists; the burden of proof is on the AKAP-positive call.
- PDE4D 267-290 and PDE2A 302-325/57-80 are correctly held below high (sensitive_only / unlikely), so this does **not** show that all PDE helices fail.

## Why this matters

This is the first real-protein confirmation of the danger-zone weakness: where an amphipathic non-AKAP also clears PSSM>=12, ML v2 provides no independent discrimination (here it actively endorsed the helix at 0.99). It is a full-protein FP in the known functional helix region — the strongest evidence class in this protocol, not a peptide-only artifact.
