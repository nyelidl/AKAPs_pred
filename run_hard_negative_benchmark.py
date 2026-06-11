#!/usr/bin/env python3
"""
run_hard_negative_benchmark.py
Runs the DEPLOYED AKAPSpred v5.1 + ML v2 on a verified hard-negative FASTA and
produces the danger-zone diagnostic. Does NOT retrain, replace, or retune anything.

Pipeline:
  1. Run akap_screen.py (CLI, exactly as deployed) with --use-ml on the FASTA.
  2. Join per-window screen output to metadata by FASTA id -> biological class.
  3. Danger zone = (amphipathic == True) AND (pssm_score >= PSSM_DZ).
  4. For each danger window, record ml_prob, proteomic_confidence, classification,
     n_negdet, is_false_positive (= proteomic_confidence in {high, very_high}),
     plus biological-context interpretation columns.
  5. Per-class summary + a written report.

Outputs:
  hard_negative_v51_screening.csv
  hard_negative_danger_zone.csv
  HARD_NEGATIVE_DIAGNOSTIC_REPORT.md

Usage:
  python3 run_hard_negative_benchmark.py \
      --fasta hard_negative_amphipathic_set.fasta \
      --meta  hard_negative_candidate_list_biological_context.csv \
      --screen /path/to/akap_screen.py     # optional; default: auto-locate
"""
import argparse, os, re, subprocess, sys, shutil, datetime
import pandas as pd

PSSM_DZ = 12.0
HIGH_TIERS = ("high", "very_high")

def locate_screen(explicit):
    if explicit and os.path.exists(explicit):
        return explicit
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, "akap_screen.py"),
              os.path.join(here, "..", "akap_screen.py"),
              "akap_screen.py"):
        if os.path.exists(c):
            return os.path.abspath(c)
    return None

def parse_region(s):
    """Best-effort 'start-end' extraction from a free-text known_region."""
    if not isinstance(s, str): return None
    m = re.search(r"(\d+)\s*[-–]\s*(\d+)", s)
    return (int(m.group(1)), int(m.group(2))) if m else None

def in_known_region(win_start, win_end, region):
    rng = parse_region(region)
    if not rng: return "unknown"
    a, b = rng
    # overlap of [win_start,win_end] with [a,b]
    return "yes" if (win_start <= b and win_end >= a) else "no"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--screen", default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    if not os.path.exists(a.fasta):
        sys.exit(f"FASTA not found: {a.fasta} (supply the verified sequences first)")
    screen = locate_screen(a.screen)
    if not screen:
        sys.exit("akap_screen.py not found; pass --screen /path/to/akap_screen.py")

    meta = pd.read_csv(a.meta)
    meta = meta[meta["id"].astype(str) != "EXAMPLE_DELETE_ME"].copy()
    meta["id"] = meta["id"].astype(str)

    # 1. run the deployed screener exactly as in production
    raw = os.path.join(a.outdir, "_raw_screen.csv")
    cmd = [sys.executable, screen, a.fasta, "--use-ml", "-o", raw]
    print("[run]", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout[-500:]); 
    if r.returncode != 0:
        print(r.stderr[-1500:]); sys.exit("akap_screen.py failed")
    if not os.path.exists(raw):
        sys.exit("no screen output produced (0 hits?). Check the FASTA.")

    screen_df = pd.read_csv(raw)

    # identify the protein-id column produced by the screener
    id_col = next((c for c in ("protein", "protein_id", "name", "id") if c in screen_df.columns), None)
    if id_col is None:
        sys.exit(f"cannot find protein-id column in screen output: {list(screen_df.columns)}")
    screen_df["_join_id"] = screen_df[id_col].astype(str).str.split().str[0].str.split("|").str[0]

    # 2. join to metadata (class + context)
    keep_meta = ["id","class","set_type","protein_or_peptide_name","cellular_location",
                 "biological_function","known_binding_partner","known_binding_surface",
                 "sequence_type","known_region","priority","intended_use","label"]
    keep_meta = [c for c in keep_meta if c in meta.columns]
    merged = screen_df.merge(meta[keep_meta], left_on="_join_id", right_on="id",
                             how="left", suffixes=("", "_meta"))
    unmatched = sorted(set(merged.loc[merged["class"].isna(), "_join_id"]))
    if unmatched:
        print(f"[warn] FASTA ids with no metadata: {unmatched}")

    merged.to_csv(os.path.join(a.outdir, "hard_negative_v51_screening.csv"), index=False)

    # 3-4. danger zone
    amphi_col = "amphipathic"
    pssm_col  = "pssm_score"
    conf_col  = "proteomic_confidence"
    for c in (amphi_col, pssm_col, conf_col):
        if c not in merged.columns:
            sys.exit(f"expected column '{c}' not in screen output")
    amphi_true = merged[amphi_col].astype(str).str.lower().isin(["true","1","yes"])
    dz = merged[amphi_true & (merged[pssm_col] >= PSSM_DZ)].copy()
    dz["is_false_positive"] = dz[conf_col].isin(HIGH_TIERS)

    ws = next((c for c in ("win_start","window_start") if c in dz.columns), None)
    we = next((c for c in ("win_end","window_end") if c in dz.columns), None)
    if ws and we and "known_region" in dz.columns:
        dz["in_known_helix_region"] = [in_known_region(s, e, reg)
                                       for s, e, reg in zip(dz[ws], dz[we], dz["known_region"])]
    else:
        dz["in_known_helix_region"] = "unknown"
    dz["is_peptide_only"] = dz.get("sequence_type", "").astype(str).eq("peptide_only")
    dz["partner_is_PKA"] = "no"   # by construction: every candidate is non_AKAP
    dz["interpretation_note"] = (
        "non-AKAP partner=" + dz.get("known_binding_partner","").astype(str)
        + "; surface=" + dz.get("known_binding_surface","").astype(str))

    dz_cols = (["_join_id","class","set_type","protein_or_peptide_name"] +
               [c for c in (ws, we, "core","window") if c and c in dz.columns] +
               [pssm_col, "ml_prob", conf_col, "classification", "n_negdet",
                "is_false_positive","in_known_helix_region","is_peptide_only",
                "partner_is_PKA","cellular_location","known_binding_partner",
                "known_binding_surface","interpretation_note"])
    dz_cols = [c for c in dz_cols if c in dz.columns]
    dz[dz_cols].to_csv(os.path.join(a.outdir, "hard_negative_danger_zone.csv"), index=False)

    # 5. per-class summary
    def cls_stats(df_all, df_dz):
        out = []
        classes = sorted(set(df_all["class"].dropna()))
        for cl in classes:
            allc = df_all[df_all["class"] == cl]
            dzc  = df_dz[df_dz["class"] == cl]
            fp   = dzc[dzc["is_false_positive"]]
            row = dict(biological_class=cl,
                       proteins=allc["_join_id"].nunique(),
                       windows_tested=len(allc),
                       danger_zone_windows=len(dzc),
                       danger_entry_rate=round(len(dzc)/max(len(allc),1), 4),
                       dz_ml_median=round(dzc["ml_prob"].median(), 3) if len(dzc) else None,
                       dz_ml_frac_ge_0p95=round((dzc["ml_prob"] >= 0.95).mean(), 3) if len(dzc) else None,
                       high_conf_FP=len(fp),
                       high_conf_FP_rate=round(len(fp)/max(len(dzc),1), 4) if len(dzc) else 0.0,
                       fp_in_known_region=int((fp["in_known_helix_region"] == "yes").sum()) if len(fp) else 0,
                       fp_peptide_only=int(fp["is_peptide_only"].sum()) if len(fp) else 0)
            out.append(row)
        return pd.DataFrame(out)

    summ = cls_stats(merged.dropna(subset=["class"]), dz.dropna(subset=["class"]))
    summ_path = os.path.join(a.outdir, "hard_negative_class_summary.csv")
    summ.to_csv(summ_path, index=False)

    # written report with the verdict logic the protocol specifies
    def verdict(r):
        if r["danger_zone_windows"] == 0:
            return "PSSM already filters this class (no danger-zone entry)."
        if r["dz_ml_frac_ge_0p95"] is not None and r["dz_ml_frac_ge_0p95"] >= 0.5:
            return "PRIORITY v3 hard negative: enters danger zone AND ML v2 endorses it."
        if r["dz_ml_median"] is not None and r["dz_ml_median"] < 0.6:
            return "ML v2 provides useful independent discrimination here."
        return "Mixed: enters danger zone, ML partial; candidate v3 hard negative."

    rep = os.path.join(a.outdir, "HARD_NEGATIVE_DIAGNOSTIC_REPORT.md")
    with open(rep, "w") as fh:
        fh.write("# Hard-negative diagnostic report — AKAPSpred v5.1 + ML v2\n\n")
        fh.write(f"_Generated {datetime.date.today().isoformat()}; model frozen, no retraining._\n\n")
        fh.write(f"Danger zone = `amphipathic == True AND pssm_score >= {PSSM_DZ}`. "
                 "Every candidate is a real non-AKAP, so any high/very_high call is a false positive.\n\n")
        fh.write("## Per-class results\n\n")
        fh.write(summ.to_markdown(index=False) if hasattr(summ, "to_markdown") else summ.to_string(index=False))
        fh.write("\n\n## Per-class verdict\n\n")
        for _, r in summ.iterrows():
            fh.write(f"- **{r['biological_class']}** — {verdict(r)}\n")
        fp_total = int(dz["is_false_positive"].sum())
        fp_region = int(((dz["is_false_positive"]) & (dz["in_known_helix_region"] == "yes")).sum())
        fp_pep = int(((dz["is_false_positive"]) & (dz["is_peptide_only"])).sum())
        fh.write("\n## Overall\n\n")
        fh.write(f"- danger-zone windows: {len(dz)}\n")
        fh.write(f"- high-confidence false positives: {fp_total}\n")
        fh.write(f"  - of which in the known functional helix region (deployment-relevant): {fp_region}\n")
        fh.write(f"  - of which peptide-only (weaker evidence, lacks protein context): {fp_pep}\n")
        fh.write("\n### Interpretation rules applied\n"
                 "- Class with no danger-zone entry -> PSSM already filters it.\n"
                 "- Danger-zone entry + low ML -> ML v2 adds real discrimination.\n"
                 "- Danger-zone entry + high ML -> priority v3 hard-negative class.\n"
                 "- FPs from full proteins in their known helix region are the strongest "
                 "evidence of a real deployment weakness; peptide-only FPs are reported separately.\n")
    print(f"\n[done] wrote:\n  hard_negative_v51_screening.csv\n  hard_negative_danger_zone.csv\n"
          f"  hard_negative_class_summary.csv\n  HARD_NEGATIVE_DIAGNOSTIC_REPORT.md")
    print("\nPer-class summary:")
    print(summ.to_string(index=False))

if __name__ == "__main__":
    main()
