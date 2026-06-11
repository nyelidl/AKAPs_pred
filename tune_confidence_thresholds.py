#!/usr/bin/env python3
"""
tune_confidence_thresholds.py
=============================
Grid-search the AKAPSpred final-confidence thresholds against the synthetic
benchmark, *without* re-running the PSSM/ML screen for every combination.

Rationale
---------
The proteomic-confidence tier is a pure function of per-hit features that the
screen already emits (classification, n_negdet, amphipathic, pssm_score,
isoform, ml_prob, background_risk).  So we run `akap_screen.py --use-ml` once
to produce a baseline prediction CSV, then re-derive the tier in memory for
each candidate threshold set.  This makes a few-thousand-cell grid search cheap.

The candidate tier function here is byte-for-byte the same logic that ships in
`akap_screen.assign_proteomic_confidence` (downgrade-not-block background
behaviour, isoform-specific PSSM thresholds, stricter ML gate under high
background risk).  Keep the two in sync.

Selection priority (from the project brief)
-------------------------------------------
  1. Keep false-positive HIGH-confidence low:
       random_background          high <= 1%
       composition_shuffled_decoy high <= 1%
       DDIP_like                  high <= 1%
  2. Improve true-positive recovery:
       positive_RII_like  high >= 60%
       positive_RI_like   high >= 70%
       positive_dual_like high >= 70%
  3. Keep negative_determinant_disrupted high as low as possible (<= 2%).

Usage
-----
  python tune_confidence_thresholds.py \
      --pred  full_pred_v5_baseline.csv \
      --labels synthetic_akap_10k_labels.csv \
      --out   threshold_tuning_summary.csv
"""
import argparse
import itertools
import pandas as pd
import numpy as np


# --------------------------------------------------------------------------- #
# Candidate confidence function (mirrors akap_screen.assign_proteomic_confidence)
# --------------------------------------------------------------------------- #
def derive_tier(row, P):
    """Return proteomic-confidence tier for one hit-row given a parameter set P."""
    classification = row["classification"]
    n_negdet       = int(row["n_negdet"]) if not pd.isna(row["n_negdet"]) else 0
    amphipathic    = bool(row["amphipathic"])
    pssm           = float(row["pssm_score"])
    iso            = row["isoform"]
    ml             = row["ml_prob"]
    bg             = row["background_risk"]
    ml = None if (ml is None or (isinstance(ml, float) and np.isnan(ml))) else float(ml)

    iso_thr = P["RII_PSSM_HIGH_THR"] if iso == "RII" else P["RI_PSSM_HIGH_THR"]

    # --- red flags ---
    if classification == "unlikely" or n_negdet >= 2:
        return "unlikely"
    if n_negdet == 1 and pssm < iso_thr:
        return "unlikely"

    bio_ok  = (classification == "AKAP" and n_negdet == 0)
    amphi_ok = amphipathic or pssm >= (iso_thr + 6)
    pssm_ok = pssm >= iso_thr

    bg_low_mod = bg in ("low", "moderate")
    bg_high    = bg in ("elevated", "high")

    ml_av = ml is not None

    if bio_ok and amphi_ok and pssm_ok:
        if bg_low_mod:
            if ml_av and ml >= P["ML_VHIGH_THR"]:
                return "very_high"
            if ml_av and ml >= P["ML_HIGH_THR"]:
                return "high"
            if not ml_av:
                return "high"          # rule-only fallback
        elif bg_high:
            # downgrade, not block: strong hit can still be HIGH (never very_high)
            pssm_floor_bg = iso_thr + P["HIGH_BG_PSSM_BONUS"]
            if ml_av and ml >= P["ML_HIGH_BG_THR"] and pssm >= pssm_floor_bg:
                return "high"
            if ml_av and ml >= P["ML_HIGH_THR"]:
                return "medium"

    if ml_av and ml >= P["ML_MEDIUM_THR"]:
        return "medium"
    return "sensitive_only"


# --------------------------------------------------------------------------- #
def evaluate(pred, labels, P):
    pred = pred.copy()
    pred["tier"] = pred.apply(lambda r: derive_tier(r, P), axis=1)
    hi_prot = set(pred.loc[pred["tier"].isin(["high", "very_high"]), "protein"].astype(str))
    vh_prot = set(pred.loc[pred["tier"] == "very_high", "protein"].astype(str))

    lab = labels[["protein_id", "category"]].copy()
    lab["protein_id"] = lab["protein_id"].astype(str)
    lab["hi"] = lab["protein_id"].isin(hi_prot)
    lab["vh"] = lab["protein_id"].isin(vh_prot)
    g = lab.groupby("category").agg(n=("protein_id", "count"),
                                    hi=("hi", "mean"),
                                    vh=("vh", "mean"))
    return {c: (g.loc[c, "hi"], g.loc[c, "vh"]) for c in g.index}


def passes_constraints(res):
    fp_ok = (res["random_background"][0]          <= 0.01 and
             res["composition_shuffled_decoy"][0] <= 0.01 and
             res["DDIP_like"][0]                   <= 0.01)
    tp_ok = (res["positive_RII_like"][0]  >= 0.60 and
             res["positive_RI_like"][0]   >= 0.70 and
             res["positive_dual_like"][0] >= 0.70)
    neg_ok = res["negative_determinant_disrupted"][0] <= 0.02
    return fp_ok, tp_ok, neg_ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", default="threshold_tuning_summary.csv")
    args = ap.parse_args()

    pred = pd.read_csv(args.pred)
    labels = pd.read_csv(args.labels)

    grid = dict(
        RII_PSSM_HIGH_THR=[10.5, 11.0, 11.5, 12.0],
        RI_PSSM_HIGH_THR=[11.5, 12.0],
        ML_HIGH_THR=[0.80, 0.85, 0.90],
        ML_HIGH_BG_THR=[0.90, 0.95, 0.97, 0.99],
        HIGH_BG_PSSM_BONUS=[0.0, 1.0, 2.0],
    )
    fixed = dict(ML_VHIGH_THR=0.90, ML_MEDIUM_THR=0.60)

    keys = list(grid.keys())
    rows = []
    for combo in itertools.product(*[grid[k] for k in keys]):
        P = dict(zip(keys, combo)); P.update(fixed)
        res = evaluate(pred, labels, P)
        fp_ok, tp_ok, neg_ok = passes_constraints(res)
        rows.append({
            **{k: P[k] for k in keys},
            "random_high":  res["random_background"][0],
            "shuffle_high": res["composition_shuffled_decoy"][0],
            "ddip_high":    res["DDIP_like"][0],
            "negdet_high":  res["negative_determinant_disrupted"][0],
            "RII_high":     res["positive_RII_like"][0],
            "RI_high":      res["positive_RI_like"][0],
            "dual_high":    res["positive_dual_like"][0],
            "RII_vhigh":    res["positive_RII_like"][1],
            "fp_ok": fp_ok, "tp_ok": tp_ok, "neg_ok": neg_ok,
            "all_ok": fp_ok and tp_ok and neg_ok,
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)

    # Selection philosophy
    # --------------------
    # The benchmark is synthetic and the ML AUC (~0.9999) is optimistic, so we
    # do NOT chase the last few points of recovery by leaning on a near-certain
    # ML gate (ML_HIGH_BG=0.99) that would not transfer to real data. We impose
    # a robustness cap (ML_HIGH_BG <= 0.95) and, among combos that still meet the
    # FP and TP targets, prefer the one with the lowest false-positive burden
    # (lowest negative_determinant_disrupted high-rate, then random FP), while
    # keeping RII recovery comfortably above target.
    ROBUST_BG_CAP = 0.95
    pool = df[df["fp_ok"] & df["tp_ok"] & (df["ML_HIGH_BG_THR"] <= ROBUST_BG_CAP)
              & (df["RII_high"] >= 0.75)].copy()
    if not len(pool):                       # relax if the cap is too tight
        pool = df[df["fp_ok"] & df["tp_ok"]].copy()
    if not len(pool):
        pool = df[df["fp_ok"]].copy()
    pool = pool.sort_values(
        by=["negdet_high", "random_high", "RII_high"],
        ascending=[True, True, False])
    best = pool.iloc[0]

    print(f"Grid cells evaluated: {len(df)}")
    print(f"Cells passing ALL constraints: {int(df['all_ok'].sum())}")
    print(f"Cells passing FP+TP: {int((df['fp_ok'] & df['tp_ok']).sum())}\n")
    print("RECOMMENDED THRESHOLD SET")
    for k in keys:
        print(f"  {k:20s} = {best[k]}")
    for k in fixed:
        print(f"  {k:20s} = {fixed[k]}")
    print("\nResulting category high-confidence rates:")
    for c in ["random_background", "composition_shuffled_decoy", "DDIP_like",
              "negative_determinant_disrupted", "positive_RII_like",
              "positive_RI_like", "positive_dual_like"]:
        key = {"random_background": "random_high",
               "composition_shuffled_decoy": "shuffle_high",
               "DDIP_like": "ddip_high",
               "negative_determinant_disrupted": "negdet_high",
               "positive_RII_like": "RII_high",
               "positive_RI_like": "RI_high",
               "positive_dual_like": "dual_high"}[c]
        print(f"  {c:34s} {best[key]*100:5.1f}%")
    print(f"\nSaved full grid to {args.out}")
    return best


if __name__ == "__main__":
    main()
