#!/usr/bin/env python3
"""
akap_ml.py
==========
Machine-learning-enhanced AKAP motif classifier.

Trains a gradient-boosted classifier on biologically motivated features
extracted from the Burgers et al. (2015) SI motif tables. The model can then
score arbitrary protein windows alongside (or instead of) the raw PSSM.

Feature set (per window)
------------------------
 1. PSSM score                         (1 feature)
 2. Hydrophobic moment  (Eisenberg, α-helix 100°)   (1)
 3. Mean hydrophobicity (Eisenberg)     (1)
 4. Helical-wheel amphipathicity index  (1)
 5. Net charge at pH 7                  (1)
 6. Fraction charged, fraction polar    (2)
 7. Helix propensity (Pace-Scholtz mean)(1)
 8. Proline count, Glycine count, Cys count (3)
 9. Amino-acid composition (20)         (20)
10. Anchor-position hydrophobicity      (n_anchors)
11. Polar-face polarity score           (1)
12. Hydrophobic moment via FFT at α-helix period (1)
                                        ───
                                        ~36-40 features

Negative set construction
-------------------------
 • Composition-shuffled SI motifs (5× per positive)
 • Random Swiss-Prot-composition windows (5× per positive)
 → balanced 1:10 positive:negative for RII, 1:10 for RI

Usage
-----
  python akap_ml.py                        # train, cross-validate, save model
  python akap_ml.py --predict seq.fasta    # predict on new sequences
  python akap_ml.py --evaluate             # detailed evaluation with plots

Requirements: scikit-learn, xgboost, pandas, numpy, openpyxl, matplotlib (optional for plots)
"""

import json
import math
import os
import random
import sys
import warnings
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, classification_report, roc_auc_score,
    precision_recall_curve, average_precision_score, roc_curve
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
import joblib

warnings.filterwarnings("ignore")

# Try optional imports
try:
    from xgboost import XGBClassifier
    _HAVE_XGB = True
except ImportError:
    _HAVE_XGB = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_PLT = True
except ImportError:
    _HAVE_PLT = False

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SI_RII = os.path.join(SCRIPT_DIR, "..", "user-data", "uploads", "bi500721a_si_001.xlsx")
SI_RI  = os.path.join(SCRIPT_DIR, "..", "user-data", "uploads", "bi500721a_si_002.xlsx")
# Fallback paths
for candidate in [
    "/mnt/user-data/uploads/bi500721a_si_001.xlsx",
    "bi500721a_si_001.xlsx",
]:
    if os.path.exists(candidate):
        SI_RII = candidate
        SI_RI  = candidate.replace("001", "002")
        break

PSSM_PATH = os.path.join(SCRIPT_DIR, "akap_pssm.json")
MODEL_PATH = os.path.join(SCRIPT_DIR, "akap_ml_model.joblib")

# ─────────────────────────────────────────────────────────────────────────────
# Physicochemical scales
# ─────────────────────────────────────────────────────────────────────────────
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

# Eisenberg consensus hydrophobicity
EISENBERG = {
    'A': 0.620, 'R':-2.530, 'N':-0.780, 'D':-0.900, 'C': 0.290,
    'Q':-0.850, 'E':-0.740, 'G': 0.480, 'H':-0.400, 'I': 1.380,
    'L': 1.060, 'K':-1.500, 'M': 0.640, 'F': 1.190, 'P': 0.120,
    'S':-0.180, 'T':-0.050, 'W': 0.810, 'Y': 0.260, 'V': 1.080,
}

PACE_SCHOLTZ = {
    'A':0.00,'L':0.21,'R':0.21,'M':0.24,'K':0.26,'Q':0.39,'E':0.40,'I':0.41,
    'W':0.49,'S':0.50,'Y':0.53,'F':0.54,'H':0.61,'V':0.61,'N':0.65,'T':0.66,
    'C':0.68,'D':0.69,'G':1.00,'P':3.16,
}

# Charge at pH 7 (simple: +1 for K/R, -1 for D/E, +0.1 for H)
CHARGE7 = {a: 0.0 for a in AA_ORDER}
CHARGE7.update({'K': 1.0, 'R': 1.0, 'H': 0.1, 'D': -1.0, 'E': -1.0})

# Swiss-Prot average amino-acid composition (for generating random negatives)
SWISSPROT_FREQ = {
    'A':0.0826,'R':0.0553,'N':0.0406,'D':0.0546,'C':0.0138,'Q':0.0393,
    'E':0.0674,'G':0.0708,'H':0.0227,'I':0.0593,'L':0.0966,'K':0.0584,
    'M':0.0242,'F':0.0386,'P':0.0470,'S':0.0657,'T':0.0534,'W':0.0108,
    'Y':0.0292,'V':0.0687,
}

POLAR_SET   = set("HRKDENQ")
CHARGED_SET = set("RKDE")
HYDRO_SET   = set("AVLIMFW")

# RII anchor positions (0-indexed within the 24-mer alignment)
RII_ANCHORS = [6, 9, 10, 13, 14, 17]
# RI anchor positions (0-indexed within the 30-mer alignment)
RI_ANCHORS  = [8, 9, 12, 13, 16, 17, 20, 21, 24, 25]

# ─────────────────────────────────────────────────────────────────────────────
# Feature extraction
# ─────────────────────────────────────────────────────────────────────────────
def hydrophobic_moment(seq, angle_deg=100.0):
    """Eisenberg hydrophobic moment for an α-helix (default 100° per residue)."""
    angle = math.radians(angle_deg)
    sin_sum = sum(EISENBERG.get(aa, 0.0) * math.sin(i * angle) for i, aa in enumerate(seq))
    cos_sum = sum(EISENBERG.get(aa, 0.0) * math.cos(i * angle) for i, aa in enumerate(seq))
    return math.sqrt(sin_sum**2 + cos_sum**2) / len(seq)


def fft_amphipathicity(seq):
    """FFT-based amphipathicity: power at the α-helix frequency (1/3.6 residues)."""
    h = np.array([EISENBERG.get(aa, 0.0) for aa in seq])
    h = h - h.mean()
    if len(h) < 4:
        return 0.0
    fft = np.fft.rfft(h)
    power = np.abs(fft)**2
    freqs = np.fft.rfftfreq(len(h))
    # α-helix period = 3.6 residues → frequency = 1/3.6 ≈ 0.278
    target = 1.0 / 3.6
    idx = np.argmin(np.abs(freqs - target))
    return float(power[idx] / power.sum()) if power.sum() > 0 else 0.0


def extract_features(seq, iso, pssm_mat=None):
    """Extract feature vector from a single window sequence."""
    n = len(seq)
    feats = {}

    # 1. PSSM score
    if pssm_mat is not None and len(pssm_mat) == n:
        feats["pssm_score"] = sum(pssm_mat[i].get(seq[i], 0.0) for i in range(n))
    else:
        feats["pssm_score"] = 0.0

    # 2. Hydrophobic moment
    feats["hydro_moment"] = hydrophobic_moment(seq)

    # 3. Mean hydrophobicity
    feats["mean_hydro"] = np.mean([EISENBERG.get(aa, 0.0) for aa in seq])

    # 4. FFT amphipathicity
    feats["fft_amphi"] = fft_amphipathicity(seq)

    # 5. Net charge
    feats["net_charge"] = sum(CHARGE7.get(aa, 0.0) for aa in seq)

    # 6. Fraction charged, fraction polar
    feats["frac_charged"] = sum(1 for aa in seq if aa in CHARGED_SET) / n
    feats["frac_polar"]   = sum(1 for aa in seq if aa in POLAR_SET) / n

    # 7. Helix propensity
    feats["helix_prop"] = np.mean([1.0 - min(PACE_SCHOLTZ.get(aa, 1.0), 1.0) for aa in seq])

    # 8. Pro, Gly, Cys counts
    feats["n_pro"] = seq.count("P")
    feats["n_gly"] = seq.count("G")
    feats["n_cys"] = seq.count("C")

    # 9. Amino-acid composition (20 features)
    counts = Counter(seq)
    for aa in AA_ORDER:
        feats[f"aa_{aa}"] = counts.get(aa, 0) / n

    # 10. Anchor-position hydrophobicity
    anchors = RII_ANCHORS if iso == "RII" else RI_ANCHORS
    anchor_h = [EISENBERG.get(seq[i], 0.0) for i in anchors if i < n]
    feats["anchor_hydro_mean"] = np.mean(anchor_h) if anchor_h else 0.0
    feats["anchor_hydro_min"]  = np.min(anchor_h) if anchor_h else 0.0

    # 11. Polar-face polarity: mean polarity at non-anchor (X) positions
    all_pos = set(range(n))
    x_pos = all_pos - set(anchors)
    x_charge = [CHARGE7.get(seq[i], 0.0) for i in sorted(x_pos) if i < n]
    feats["polar_face_charge"] = np.mean(x_charge) if x_charge else 0.0

    # 12. Hydrophobicity contrast (anchor mean - non-anchor mean)
    non_anchor_h = [EISENBERG.get(seq[i], 0.0) for i in sorted(x_pos) if i < n]
    feats["hydro_contrast"] = feats["anchor_hydro_mean"] - (np.mean(non_anchor_h) if non_anchor_h else 0.0)

    return feats


def feats_to_array(feat_dict_list):
    """Convert list of feature dicts to (X, feature_names)."""
    names = list(feat_dict_list[0].keys())
    X = np.array([[d[k] for k in names] for d in feat_dict_list])
    return X, names


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────
def load_si_seqs(xlsx_path, sheet_name, n_cols):
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=2).dropna(subset=["Number"])
    cols = sorted(
        [c for c in df.columns if str(c).replace(".0","").isdigit()],
        key=lambda x: int(float(x))
    )[:n_cols]
    seqs = []
    for _, row in df.iterrows():
        seq = "".join(str(row[c]).strip() for c in cols)
        seqs.append(seq)
    return seqs


def generate_negatives(positives, n_per_pos=5, seed=42):
    """Generate negatives: half shuffled, half random composition."""
    rng = random.Random(seed)
    aa_pool = list(SWISSPROT_FREQ.keys())
    aa_weights = [SWISSPROT_FREQ[a] for a in aa_pool]
    negs = []
    wlen = len(positives[0])
    # Shuffled
    for seq in positives:
        for _ in range(n_per_pos // 2 + 1):
            s = list(seq)
            rng.shuffle(s)
            negs.append("".join(s))
    # Random composition
    n_random = len(positives) * (n_per_pos - n_per_pos // 2 - 1 + 1)
    for _ in range(max(n_random, len(positives) * 3)):
        s = rng.choices(aa_pool, weights=aa_weights, k=wlen)
        negs.append("".join(s))
    return negs[:len(positives) * n_per_pos]


# ─────────────────────────────────────────────────────────────────────────────
# Training
# ─────────────────────────────────────────────────────────────────────────────
def train_and_evaluate(iso, pos_seqs, pssm_mat, n_neg_ratio=10):
    """Train a classifier for one isoform, cross-validate, return model + metrics."""
    print(f"\n{'═'*60}")
    print(f"  Training {iso} classifier  ({len(pos_seqs)} positives)")
    print(f"{'═'*60}")

    # Generate negatives
    neg_seqs = generate_negatives(pos_seqs, n_per_pos=n_neg_ratio)
    print(f"  Generated {len(neg_seqs)} negatives ({n_neg_ratio}× ratio)")

    # Extract features
    print("  Extracting features...")
    pos_feats = [extract_features(s, iso, pssm_mat) for s in pos_seqs]
    neg_feats = [extract_features(s, iso, pssm_mat) for s in neg_seqs]

    all_feats = pos_feats + neg_feats
    y = np.array([1]*len(pos_feats) + [0]*len(neg_feats))
    X, feat_names = feats_to_array(all_feats)

    # Scale
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Models to try
    models = {
        "GradientBoosting": GradientBoostingClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.1,
            subsample=0.8, random_state=42),
        "RandomForest": RandomForestClassifier(
            n_estimators=300, max_depth=6, class_weight="balanced",
            random_state=42),
    }
    if _HAVE_XGB:
        scale_pos = len(neg_feats) / len(pos_feats)
        models["XGBoost"] = XGBClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.1,
            scale_pos_weight=scale_pos, random_state=42,
            use_label_encoder=False, eval_metric="logloss", verbosity=0)

    # Cross-validate each
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    best_name, best_auc, best_model = None, 0.0, None

    for name, model in models.items():
        proba = cross_val_predict(model, X_scaled, y, cv=cv, method="predict_proba")[:, 1]
        preds = (proba >= 0.5).astype(int)
        auc = roc_auc_score(y, proba)
        ap  = average_precision_score(y, proba)
        acc = accuracy_score(y, preds)

        # Recall on positives only
        pos_mask = y == 1
        recall_pos = (preds[pos_mask] == 1).mean()

        print(f"\n  {name}:")
        print(f"    ROC-AUC  = {auc:.4f}")
        print(f"    Avg-Prec = {ap:.4f}")
        print(f"    Accuracy = {acc:.4f}")
        print(f"    Recall (positives) = {recall_pos:.4f} ({int(recall_pos*pos_mask.sum())}/{pos_mask.sum()})")

        if auc > best_auc:
            best_auc, best_name, best_model = auc, name, model
            best_proba = proba

    # Retrain best on all data
    print(f"\n  Best model: {best_name} (AUC={best_auc:.4f})")
    best_model.fit(X_scaled, y)

    # Feature importances
    if hasattr(best_model, "feature_importances_"):
        imp = best_model.feature_importances_
        idx = np.argsort(imp)[::-1]
        print(f"\n  Top 15 features:")
        for rank, i in enumerate(idx[:15], 1):
            print(f"    {rank:2d}. {feat_names[i]:25s}  importance={imp[i]:.4f}")

    # Detailed recall analysis
    print(f"\n  Recall at different probability thresholds:")
    for thr in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        pos_proba = best_proba[y == 1]
        neg_proba = best_proba[y == 0]
        rec = (pos_proba >= thr).mean()
        fpr = (neg_proba >= thr).mean()
        print(f"    P>={thr:.1f}:  recall={rec:.3f}  FPR={fpr:.4f}")

    # Missed positives at P>=0.5
    missed_idx = np.where((y == 1) & (best_proba < 0.5))[0]
    if len(missed_idx) > 0:
        print(f"\n  {len(missed_idx)} positive(s) missed at P>=0.5:")
        for mi in missed_idx[:10]:
            print(f"    seq: {pos_seqs[mi]}  P={best_proba[mi]:.3f}  pssm={pos_feats[mi]['pssm_score']:.1f}")

    return dict(
        model=best_model, scaler=scaler, feat_names=feat_names,
        iso=iso, best_name=best_name, auc=best_auc,
        y=y, proba=best_proba
    )


def save_models(results_rii, results_ri, path=MODEL_PATH):
    """Save trained models + scalers to disk."""
    bundle = {}
    for res in [results_rii, results_ri]:
        if res is None:
            continue
        iso = res["iso"]
        bundle[iso] = {
            "model": res["model"],
            "scaler": res["scaler"],
            "feat_names": res["feat_names"],
            "best_name": res["best_name"],
            "auc": res["auc"],
        }
    joblib.dump(bundle, path)
    print(f"\n  Models saved to {path} ({os.path.getsize(path)//1024} KB)")


# ─────────────────────────────────────────────────────────────────────────────
# Prediction on new sequences
# ─────────────────────────────────────────────────────────────────────────────
def predict_window(seq, iso, bundle, pssm_mat):
    """Score a single window. Returns probability of being a real AKAP motif."""
    feats = extract_features(seq, iso, pssm_mat)
    X = np.array([[feats[k] for k in bundle["feat_names"]]])
    X_scaled = bundle["scaler"].transform(X)
    prob = bundle["model"].predict_proba(X_scaled)[0, 1]
    return float(prob), feats


def scan_protein_ml(seq, iso, bundle, pssm, win_len, threshold=0.5):
    """Slide a window over a protein and return ML-scored hits."""
    pssm_mat = pssm[iso]["pssm"]
    hits = []
    if len(seq) < win_len:
        return hits
    for i in range(len(seq) - win_len + 1):
        win = seq[i:i+win_len]
        prob, feats = predict_window(win, iso, bundle, pssm_mat)
        if prob >= threshold:
            hits.append(dict(
                start=i+1, end=i+win_len, window=win,
                ml_prob=round(prob, 4),
                pssm_score=round(feats["pssm_score"], 2),
                hydro_moment=round(feats["hydro_moment"], 4),
                fft_amphi=round(feats["fft_amphi"], 4),
            ))
    # Non-maximum suppression
    hits.sort(key=lambda h: -h["ml_prob"])
    chosen = []
    for h in hits:
        if all(abs(h["start"] - c["start"]) >= win_len for c in chosen):
            chosen.append(h)
    return sorted(chosen, key=lambda h: h["start"])


# ─────────────────────────────────────────────────────────────────────────────
# Plotting (if matplotlib available)
# ─────────────────────────────────────────────────────────────────────────────
def plot_results(results_rii, results_ri, out_dir=SCRIPT_DIR):
    if not _HAVE_PLT:
        print("\n  [skip] matplotlib not available — no plots generated.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("AKAP ML Classifier Validation", fontsize=14, fontweight="bold")

    for col, (res, label) in enumerate([(results_rii, "PKA-RIIα"), (results_ri, "PKA-RIα")]):
        if res is None:
            continue
        y, proba = res["y"], res["proba"]

        # ROC
        ax = axes[0, col]
        fpr, tpr, _ = roc_curve(y, proba)
        ax.plot(fpr, tpr, "b-", lw=2, label=f'AUC={res["auc"]:.3f}')
        ax.plot([0,1], [0,1], "k--", alpha=0.3)
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_title(f"{label} ROC ({res['best_name']})")
        ax.legend()

        # Score distribution
        ax = axes[1, col]
        pos_p = proba[y == 1]
        neg_p = proba[y == 0]
        ax.hist(neg_p, bins=50, alpha=0.5, color="gray", label="Negative", density=True)
        ax.hist(pos_p, bins=50, alpha=0.7, color="steelblue", label="Positive", density=True)
        ax.axvline(0.5, color="red", ls="--", label="Threshold=0.5")
        ax.set_xlabel("ML Probability")
        ax.set_ylabel("Density")
        ax.set_title(f"{label} Score Distributions")
        ax.legend()

    plt.tight_layout()
    plot_path = os.path.join(out_dir, "akap_ml_validation.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"\n  Plot saved to {plot_path}")
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    import argparse
    ap = argparse.ArgumentParser(description="ML-enhanced AKAP motif classifier")
    ap.add_argument("--predict", help="FASTA file to predict on (uses saved model)")
    ap.add_argument("--ml-thr", type=float, default=0.5, help="ML probability threshold (default 0.5)")
    ap.add_argument("--evaluate", action="store_true", help="run full evaluation with plots")
    ap.add_argument("-o", "--out", help="output CSV for predictions")
    args = ap.parse_args()

    # ── Prediction mode ──
    if args.predict:
        if not os.path.exists(MODEL_PATH):
            sys.exit(f"ERROR: no saved model at {MODEL_PATH}. Run without --predict first to train.")
        bundle = joblib.load(MODEL_PATH)
        with open(PSSM_PATH) as f:
            pssm = json.load(f)

        sys.path.insert(0, SCRIPT_DIR)
        import akap_screen as A
        import csv

        rows = []
        for pid, seq in A.read_fasta(args.predict):
            for iso, wlen in [("RII", 24), ("RI", 30)]:
                if iso not in bundle:
                    continue
                hits = scan_protein_ml(seq, iso, bundle[iso], pssm, wlen, args.ml_thr)
                for h in hits:
                    h["protein"] = pid
                    h["isoform"] = iso
                    rows.append(h)

        rows.sort(key=lambda r: -r["ml_prob"])
        fields = ["protein","isoform","start","end","window","ml_prob","pssm_score",
                  "hydro_moment","fft_amphi"]
        out = open(args.out, "w", newline="") if args.out else sys.stdout
        w = csv.DictWriter(out, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})
        if args.out:
            out.close()
            print(f"  {len(rows)} hit(s) -> {args.out}")
        return

    # ── Training mode ──
    print("="*60)
    print("  AKAP ML Classifier — Training & Validation")
    print("="*60)

    # Load PSSM
    with open(PSSM_PATH) as f:
        pssm = json.load(f)

    # Load SI data
    rii_seqs = load_si_seqs(SI_RII, "PKA-RIIa Motifs", 24)
    ri_seqs  = load_si_seqs(SI_RI,  "PKA-RIa Motifs", 30)
    print(f"\n  Loaded {len(rii_seqs)} RII and {len(ri_seqs)} RI positive motifs")

    # Train
    res_rii = train_and_evaluate("RII", rii_seqs, pssm["RII"]["pssm"], n_neg_ratio=10)
    res_ri  = train_and_evaluate("RI",  ri_seqs,  pssm["RI"]["pssm"],  n_neg_ratio=10)

    # Save
    save_models(res_rii, res_ri)

    # Plot
    if args.evaluate or _HAVE_PLT:
        plot_results(res_rii, res_ri)

    # ── Quick self-test ──
    print(f"\n{'═'*60}")
    print("  Quick self-test on known AKAP peptides")
    print(f"{'═'*60}")
    bundle = joblib.load(MODEL_PATH)
    test_cases = {
        "AKAP7g(RII)":       ("AELVRLSKRLVENAVLKAVQQYLE", "RII", 24),
        "AKAP10(RII)":       ("EAQEELAWKIAKMIVSDIMQQAQY", "RII", 24),
        "vlAKAP(RII)":       ("CLLEDKARELVNEIIYVAQEKLRN", "RII", 24),
        "smAKAP(RI)":        ("GTNTVILEYAHRLSQDILCDALQQWACNNI", "RI", 30),
        "Ezrin_Nterm(RI)":   ("RAKFYPEDVAEELIQDITQKLFFLQVKEGI", "RI", 30),
        "shuffled_neg":      ("YQADQMVAKEIESLKIIMVDAQWQ", "RII", 24),
    }
    for name, (seq, iso, wlen) in test_cases.items():
        if iso in bundle:
            prob, feats = predict_window(seq, iso, bundle[iso], pssm[iso]["pssm"])
            print(f"  {name:22s}  P={prob:.3f}  pssm={feats['pssm_score']:.1f}  "
                  f"μH={feats['hydro_moment']:.3f}  {'✓ HIT' if prob>=0.5 else '✗ miss'}")

    print(f"\n{'═'*60}")
    print("  Done! Use --predict <fasta> to score new proteins.")
    print(f"{'═'*60}\n")


if __name__ == "__main__":
    main()
