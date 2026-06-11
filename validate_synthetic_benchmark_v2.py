#!/usr/bin/env python3
"""Validate AKAPSpred synthetic benchmark predictions.

Reads labels + predictions and reports protein-level recovery/false-positive
rates for sensitive, ML, high-confidence, and very-high-confidence filters.
Also generates simple matplotlib plots.
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def _truth_category(cat: str) -> int:
    return 1 if str(cat).startswith("positive_") else 0


# --------------------------------------------------------------------------- #
# Window-level on-target / off-target analysis
# --------------------------------------------------------------------------- #
def _add_window_overlap(mpred: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    """
    For every predicted hit, compute overlap with the protein's embedded motif.

    Label coordinates `motif_start`/`motif_end` are treated as 0-based half-open
    (start inclusive, end exclusive) as written by the synthetic generator, i.e.
    the motif occupies sequence positions [motif_start+1 .. motif_end] in the
    1-based coordinates that akap_screen emits for win_start/win_end.

    Adds columns:
      overlaps_embedded_motif, overlap_length, overlap_fraction,
      on_target_hit, off_target_hit
    """
    lab = labels[["protein_id", "motif_start", "motif_end", "embedded_window"]].copy()
    lab["protein_id"] = lab["protein_id"].astype(str)
    m = mpred.merge(lab, on="protein_id", how="left", suffixes=("", "_lab"))

    has_motif = m["motif_start"].notna() & m["motif_end"].notna()

    # 1-based inclusive motif span
    ms = pd.to_numeric(m["motif_start"], errors="coerce") + 1   # ->1-based inclusive start
    me = pd.to_numeric(m["motif_end"], errors="coerce")         # already inclusive end
    ws = pd.to_numeric(m["win_start"], errors="coerce")
    we = pd.to_numeric(m["win_end"], errors="coerce")

    inter = (np.minimum(we, me) - np.maximum(ws, ms) + 1).clip(lower=0)
    inter = inter.where(has_motif, other=np.nan)

    motif_len = (me - ms + 1).where(has_motif)
    win_len   = (we - ws + 1)
    denom = pd.concat([motif_len, win_len], axis=1).min(axis=1)

    m["overlap_length"] = inter
    m["overlap_fraction"] = (inter / denom).where(has_motif)
    m["overlaps_embedded_motif"] = (inter > 0).where(has_motif, other=False).fillna(False).astype(bool)

    # A hit is "high confidence" if it is in high/very_high tier
    conf = m.get("proteomic_confidence", pd.Series([""] * len(m)))
    m["_hit_high"] = conf.isin(["high", "very_high"])

    # on/off target only meaningful where a motif exists
    m["on_target_hit"]  = (m["overlaps_embedded_motif"] & has_motif).fillna(False).astype(bool)
    m["off_target_hit"] = ((~m["overlaps_embedded_motif"]) & has_motif).fillna(False).astype(bool)
    return m


def _window_level_summary(mwin: pd.DataFrame) -> pd.DataFrame:
    """
    Per-category window-level rates, distinguishing on-target vs off-target hits
    and protein-level high-confidence attribution.
    """
    rows = []
    for cat, g in mwin.groupby("category"):
        n_hits = len(g)
        on = g[g["on_target_hit"]]
        off = g[g["off_target_hit"]]
        # protein-level: among proteins of this category that have any hit,
        # how many have an on-target hit / on-target HIGH-confidence hit / etc.
        prot = g.groupby("protein")
        any_on  = prot["on_target_hit"].any()
        any_off = prot["off_target_hit"].any()
        on_high  = g[g["on_target_hit"]  & g["_hit_high"]].groupby("protein").size().reindex(any_on.index, fill_value=0) > 0
        off_high = g[g["off_target_hit"] & g["_hit_high"]].groupby("protein").size().reindex(any_on.index, fill_value=0) > 0
        n_prot_with_hit = g["protein"].nunique()
        rows.append({
            "category": cat,
            "n_proteins_with_hit": n_prot_with_hit,
            "n_windows": n_hits,
            "n_on_target_windows": int(g["on_target_hit"].sum()),
            "n_off_target_windows": int(g["off_target_hit"].sum()),
            "on_target_any_hit_rate":  float(any_on.mean()) if len(any_on) else np.nan,
            "off_target_any_hit_rate": float(any_off.mean()) if len(any_off) else np.nan,
            "on_target_high_confidence_rate":  float(on_high.mean()) if len(on_high) else np.nan,
            "off_target_high_confidence_rate": float(off_high.mean()) if len(off_high) else np.nan,
        })
    return pd.DataFrame(rows)


def summarize(labels: pd.DataFrame, pred: pd.DataFrame):
    labels = labels.copy()
    labels["is_positive"] = labels["category"].map(_truth_category)

    # Normalize prediction column names
    if "protein" not in pred.columns:
        raise ValueError("Prediction CSV must contain a 'protein' column")

    pred = pred.copy()
    if "passes_sensitive_filter" not in pred.columns:
        pred["passes_sensitive_filter"] = True
    if "passes_ml_filter" not in pred.columns:
        pred["passes_ml_filter"] = False
    if "passes_proteomic_filter" not in pred.columns:
        pred["passes_proteomic_filter"] = pred.get("proteomic_confidence", "") .isin(["high", "very_high"])
    if "passes_very_high_filter" not in pred.columns:
        pred["passes_very_high_filter"] = pred.get("proteomic_confidence", "") .eq("very_high")
    if "ml_prob" not in pred.columns:
        pred["ml_prob"] = np.nan

    protein_ids = labels["protein_id"].astype(str)
    any_hit = set(pred["protein"].astype(str).unique())
    sensitive = set(pred.loc[pred["passes_sensitive_filter"].astype(bool), "protein"].astype(str))
    ml_pass = set(pred.loc[pred["passes_ml_filter"].astype(bool), "protein"].astype(str))
    high = set(pred.loc[pred["passes_proteomic_filter"].astype(bool), "protein"].astype(str))
    very_high = set(pred.loc[pred["passes_very_high_filter"].astype(bool), "protein"].astype(str))

    prot = labels[["protein_id", "category", "is_positive"]].copy()
    prot["protein_id"] = prot["protein_id"].astype(str)
    prot["any_hit"] = prot["protein_id"].isin(any_hit)
    prot["sensitive_pass"] = prot["protein_id"].isin(sensitive)
    prot["ml_pass"] = prot["protein_id"].isin(ml_pass)
    prot["high_confidence"] = prot["protein_id"].isin(high)
    prot["very_high_confidence"] = prot["protein_id"].isin(very_high)

    cat = prot.groupby("category").agg(
        n=("protein_id", "count"),
        any_hit=("any_hit", "sum"),
        sensitive_pass=("sensitive_pass", "sum"),
        ml_pass=("ml_pass", "sum"),
        high_confidence=("high_confidence", "sum"),
        very_high_confidence=("very_high_confidence", "sum"),
    ).reset_index()
    for col in ["any_hit", "sensitive_pass", "ml_pass", "high_confidence", "very_high_confidence"]:
        cat[col + "_rate"] = cat[col] / cat["n"]

    # Overall positive recovery and decoy/background FPRs
    positives = prot[prot["is_positive"] == 1]
    negatives = prot[prot["is_positive"] == 0]
    metrics = {
        "n_proteins": len(prot),
        "n_predictions": len(pred),
        "positive_recovery_any_hit": float(positives["any_hit"].mean()) if len(positives) else np.nan,
        "positive_recovery_high": float(positives["high_confidence"].mean()) if len(positives) else np.nan,
        "positive_recovery_very_high": float(positives["very_high_confidence"].mean()) if len(positives) else np.nan,
        "negative_fpr_any_hit": float(negatives["any_hit"].mean()) if len(negatives) else np.nan,
        "negative_fpr_ml_pass": float(negatives["ml_pass"].mean()) if len(negatives) else np.nan,
        "negative_fpr_high": float(negatives["high_confidence"].mean()) if len(negatives) else np.nan,
        "negative_fpr_very_high": float(negatives["very_high_confidence"].mean()) if len(negatives) else np.nan,
    }
    return prot, cat, metrics, pred.merge(labels[["protein_id", "category", "is_positive"]], left_on="protein", right_on="protein_id", how="left")


def make_plots(cat: pd.DataFrame, mpred: pd.DataFrame, outdir: Path):
    import matplotlib.pyplot as plt
    outdir.mkdir(parents=True, exist_ok=True)

    plot_df = cat.sort_values("category")
    x = np.arange(len(plot_df))
    width = 0.20
    plt.figure(figsize=(14, 5))
    for i, col in enumerate(["any_hit_rate", "ml_pass_rate", "high_confidence_rate", "very_high_confidence_rate"]):
        plt.bar(x + (i - 1.5) * width, plot_df[col], width, label=col.replace("_", " "))
    plt.xticks(x, plot_df["category"], rotation=45, ha="right")
    plt.ylabel("Protein-level rate")
    plt.ylim(0, 1.05)
    plt.title("AKAPSpred synthetic benchmark: PSSM vs ML vs final filters")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "filter_comparison_by_category.png", dpi=200)
    plt.close()

    if "ml_prob" in mpred.columns and mpred["ml_prob"].notna().any():
        plt.figure(figsize=(12, 5))
        cats = [c for c in sorted(mpred["category"].dropna().unique())]
        data = [mpred.loc[mpred["category"] == c, "ml_prob"].dropna().values for c in cats]
        if any(len(d) for d in data):
            plt.boxplot(data, labels=cats, showfliers=False)
            plt.xticks(rotation=45, ha="right")
            plt.ylabel("ML probability")
            plt.title("ML probability distribution by synthetic category")
            plt.tight_layout()
            plt.savefig(outdir / "ml_probability_by_category.png", dpi=200)
        plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True, help="synthetic_akap_10k_labels.csv")
    ap.add_argument("--pred", required=True, help="AKAPSpred prediction CSV")
    ap.add_argument("--outdir", default="validation_results_v2")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    labels = pd.read_csv(args.labels)
    pred = pd.read_csv(args.pred)

    prot, cat, metrics, mpred = summarize(labels, pred)

    # ── Window-level on-target / off-target analysis ──
    mwin = _add_window_overlap(mpred, labels)
    win_summary = _window_level_summary(mwin)
    # fold on/off-target high-confidence rates into the category summary
    cat = cat.merge(
        win_summary[["category", "on_target_any_hit_rate", "on_target_high_confidence_rate",
                     "off_target_any_hit_rate", "off_target_high_confidence_rate"]],
        on="category", how="left")

    prot.to_csv(outdir / "protein_level_summary.csv", index=False)
    cat.to_csv(outdir / "category_summary.csv", index=False)
    pd.DataFrame([metrics]).to_csv(outdir / "overall_metrics.csv", index=False)
    mwin.to_csv(outdir / "predictions_with_window_overlap.csv", index=False)
    win_summary.to_csv(outdir / "window_level_summary.csv", index=False)
    make_plots(cat, mpred, outdir)

    print("\nCategory summary (protein-level + on/off-target high-confidence):")
    show_cols = ["category", "n", "any_hit_rate", "ml_pass_rate",
                 "high_confidence_rate", "very_high_confidence_rate",
                 "on_target_high_confidence_rate", "off_target_high_confidence_rate"]
    show_cols = [c for c in show_cols if c in cat.columns]
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(cat[show_cols].to_string(index=False))
    print("\nWindow-level summary:")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(win_summary.to_string(index=False))
    print("\nOverall metrics:")
    for k, v in metrics.items():
        print(f"  {k}: {v}")
    print(f"\nSaved validation outputs to: {outdir}")


if __name__ == "__main__":
    main()
