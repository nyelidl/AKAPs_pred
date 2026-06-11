#!/usr/bin/env python3
"""
train_akap_ml_v2.py
===================
Train an improved ML classifier for AKAP motif discrimination.

Key improvements over v1:
  1. Hard-negative mining: PSSM-positive false positives from shuffled/random sequences
  2. Multi-model comparison: LR, RF, GBT, XGBoost
  3. Selects model with lowest FPR on hard negatives while retaining positive recall
  4. Saves: akap_ml_model_v2.joblib + akap_ml_model_v2_metadata.json

Usage:
  python train_akap_ml_v2.py                          # default paths
  python train_akap_ml_v2.py --si-rii path/001.xlsx   # custom SI paths
"""

import argparse, json, math, os, random, sys, warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, classification_report
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
import joblib

warnings.filterwarnings("ignore")

try:
    from xgboost import XGBClassifier
    _HAVE_XGB = True
except ImportError:
    _HAVE_XGB = False

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from akap_ml import extract_features, AA_ORDER, SWISSPROT_FREQ


def load_si_seqs(xlsx_path, sheet_name, n_cols):
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=2).dropna(subset=["Number"])
    cols = sorted([c for c in df.columns if str(c).replace(".0","").isdigit()],
                  key=lambda x: int(float(x)))[:n_cols]
    return ["".join(str(r[c]).strip() for c in cols) for _, r in df.iterrows()]


def generate_easy_negatives(positives, n_per_pos=3, seed=42):
    """Shuffled + random-composition negatives."""
    rng = random.Random(seed)
    aa_pool = list(SWISSPROT_FREQ.keys())
    aa_w = [SWISSPROT_FREQ[a] for a in aa_pool]
    negs = []
    wlen = len(positives[0])
    for seq in positives:
        for _ in range(n_per_pos // 2 + 1):
            s = list(seq); rng.shuffle(s); negs.append("".join(s))
    n_rand = len(positives) * (n_per_pos - n_per_pos // 2)
    for _ in range(n_rand):
        negs.append("".join(rng.choices(aa_pool, weights=aa_w, k=wlen)))
    return negs[:len(positives) * n_per_pos]


def generate_hard_negatives(positives, pssm, iso, n_target=500, seed=42):
    """
    Hard negatives: random sequences that PASS the PSSM threshold.
    These are the ones that cause false positives in proteomic screening.
    """
    from akap_screen import scan_isoform, DEFAULT_RII_THR, DEFAULT_RI_THR
    rng = random.Random(seed)
    aa_pool = list(SWISSPROT_FREQ.keys())
    aa_w = [SWISSPROT_FREQ[a] for a in aa_pool]
    wlen = 24 if iso == "RII" else 30
    thr = DEFAULT_RII_THR if iso == "RII" else DEFAULT_RI_THR
    hard_negs = []
    attempts = 0
    max_attempts = n_target * 200  # generous limit
    while len(hard_negs) < n_target and attempts < max_attempts:
        # Generate a longer random protein (200-800 aa)
        prot_len = rng.randint(200, 800)
        prot = "".join(rng.choices(aa_pool, weights=aa_w, k=prot_len))
        hits = scan_isoform(prot, iso, pssm, thr)
        for h in hits:
            hard_negs.append(h["window"])
            if len(hard_negs) >= n_target:
                break
        attempts += 1
    print(f"  Generated {len(hard_negs)} hard negatives for {iso} ({attempts} attempts)")
    return hard_negs


def generate_negdet_negatives(positives, n_target=200, seed=42):
    """Create negatives by inserting D/E at hydrophobic anchor positions."""
    rng = random.Random(seed)
    negs = []
    for seq in positives:
        s = list(seq)
        # Replace 2-3 hydrophobic-face positions with D or E
        positions = [6, 9, 10, 13, 14, 17] if len(seq) == 24 else [8, 12, 16, 20, 24]
        for pos in rng.sample(positions, min(3, len(positions))):
            if pos < len(s):
                s[pos] = rng.choice("DE")
        negs.append("".join(s))
        if len(negs) >= n_target:
            break
    return negs


def main():
    ap = argparse.ArgumentParser(description="Train AKAP ML v2 with hard-negative mining")
    ap.add_argument("--si-rii", default=None, help="path to SI xlsx (RII motifs)")
    ap.add_argument("--si-ri", default=None, help="path to SI xlsx (RI motifs)")
    ap.add_argument("--out", default=os.path.join(SCRIPT_DIR, "akap_ml_model_v2.joblib"))
    args = ap.parse_args()

    # Find SI files
    si_rii = args.si_rii
    si_ri = args.si_ri
    for d in [SCRIPT_DIR, "/mnt/user-data/uploads", ".", "/content"]:
        if si_rii is None and os.path.exists(os.path.join(d, "bi500721a_si_001.xlsx")):
            si_rii = os.path.join(d, "bi500721a_si_001.xlsx")
            si_ri = os.path.join(d, "bi500721a_si_002.xlsx")
    if si_rii is None:
        sys.exit("ERROR: cannot find SI xlsx files. Use --si-rii and --si-ri.")

    import akap_screen as A
    pssm_data = A.load_pssm()

    print("=" * 60)
    print("  AKAP ML v2 Training — Hard-Negative Mining")
    print("=" * 60)

    results = {}
    for iso, si_path, sheet, ncol in [
        ("RII", si_rii, "PKA-RIIa Motifs", 24),
        ("RI", si_ri, "PKA-RIa Motifs", 30),
    ]:
        print(f"\n{'─'*60}")
        print(f"  Training {iso} classifier")
        print(f"{'─'*60}")

        pos_seqs = load_si_seqs(si_path, sheet, ncol)
        print(f"  Positives: {len(pos_seqs)}")

        # Generate negatives
        easy_negs = generate_easy_negatives(pos_seqs, n_per_pos=3)
        print(f"  Easy negatives: {len(easy_negs)}")

        hard_negs = generate_hard_negatives(pos_seqs, pssm_data, iso, n_target=min(500 if iso=="RII" else 100, 500))
        negdet_negs = generate_negdet_negatives(pos_seqs, n_target=200)
        print(f"  Neg-det disrupted negatives: {len(negdet_negs)}")

        all_negs = easy_negs + hard_negs + negdet_negs
        print(f"  Total negatives: {len(all_negs)}")

        # Extract features
        pssm_mat = pssm_data[iso]["pssm"]
        print("  Extracting features...")
        pos_feats = [extract_features(s, iso, pssm_mat) for s in pos_seqs]
        neg_feats = [extract_features(s, iso, pssm_mat) for s in all_negs]

        feat_names = list(pos_feats[0].keys())
        X = np.array([[d[k] for k in feat_names] for d in pos_feats + neg_feats])
        y = np.array([1]*len(pos_feats) + [0]*len(neg_feats))

        # Mark hard negatives for evaluation
        n_pos = len(pos_feats)
        n_easy = len(easy_negs)
        n_hard = len(hard_negs)
        n_negdet = len(negdet_negs)

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # ── Multi-model comparison ──
        models = {
            "LogisticRegression": LogisticRegression(max_iter=1000, C=1.0, random_state=42),
            "RandomForest": RandomForestClassifier(
                n_estimators=300, max_depth=6, class_weight="balanced", random_state=42),
            "GradientBoosting": GradientBoostingClassifier(
                n_estimators=200, max_depth=4, learning_rate=0.1, subsample=0.8, random_state=42),
        }
        if _HAVE_XGB:
            scale_pos = len(all_negs) / len(pos_seqs)
            models["XGBoost"] = XGBClassifier(
                n_estimators=200, max_depth=4, learning_rate=0.1,
                scale_pos_weight=scale_pos, random_state=42,
                use_label_encoder=False, eval_metric="logloss", verbosity=0)

        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        best_name, best_score, best_model = None, -1, None

        for name, model in models.items():
            proba = cross_val_predict(model, X_scaled, y, cv=cv, method="predict_proba")[:, 1]

            auc = roc_auc_score(y, proba)
            ap_score = average_precision_score(y, proba)

            # Key metric: FPR on hard negatives at threshold 0.80
            hard_start = n_pos + n_easy
            hard_end = hard_start + n_hard
            hard_proba = proba[hard_start:hard_end]
            hard_fpr_80 = (hard_proba >= 0.80).mean() if len(hard_proba) > 0 else 0

            # Positive recall at 0.80
            pos_proba = proba[:n_pos]
            pos_recall_80 = (pos_proba >= 0.80).mean()

            # Composite score: maximize recall while minimizing hard-neg FPR
            composite = pos_recall_80 - 2 * hard_fpr_80 + 0.1 * auc

            print(f"\n  {name}:")
            print(f"    ROC-AUC={auc:.4f}  PR-AUC={ap_score:.4f}")
            print(f"    Pos recall @0.80: {pos_recall_80:.3f}")
            print(f"    Hard-neg FPR @0.80: {hard_fpr_80:.3f}")
            print(f"    Composite score: {composite:.4f}")

            if composite > best_score:
                best_score, best_name, best_model = composite, name, model

        print(f"\n  → Best model: {best_name} (composite={best_score:.4f})")

        # Retrain on full data
        best_model.fit(X_scaled, y)

        # Feature importances
        if hasattr(best_model, "feature_importances_"):
            imp = best_model.feature_importances_
            top = np.argsort(imp)[::-1][:10]
            print(f"  Top features:")
            for rank, i in enumerate(top[:5], 1):
                print(f"    {rank}. {feat_names[i]:25s} {imp[i]:.4f}")

        results[iso] = {
            "model": best_model,
            "scaler": scaler,
            "feat_names": feat_names,
            "best_name": best_name,
            "auc": roc_auc_score(y, best_model.predict_proba(X_scaled)[:, 1]),
            "n_pos": n_pos,
            "n_neg_total": len(all_negs),
            "n_hard_neg": n_hard,
        }

    # Save
    joblib.dump(results, args.out)
    meta = {iso: {k: v for k, v in d.items() if k != "model" and k != "scaler"}
            for iso, d in results.items()}
    meta_path = args.out.replace(".joblib", "_metadata.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2, default=str)

    print(f"\n{'═'*60}")
    print(f"  Saved: {args.out}")
    print(f"  Saved: {meta_path}")
    print(f"{'═'*60}")


if __name__ == "__main__":
    main()
