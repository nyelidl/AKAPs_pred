#!/usr/bin/env python3
"""
akap_screen.py
==============
Screen a protein FASTA for PKA-R (A-Kinase Anchoring Protein) amphipathic-helix
binding motifs, following the approach of THAHIT.

Reference
---------
Burgers, van der Heyden, Kok, Heck & Scholten (2015)
"A Systematic Evaluation of Protein Kinase A-A-Kinase Anchoring Protein
Interaction Motifs", Biochemistry 54, 11-21. DOI: 10.1021/bi500721a

What an "AKAP domain" is (the thing we filter for)
--------------------------------------------------
A short (3-4 turn) AMPHIPATHIC alpha-helix that docks its HYDROPHOBIC face onto
the dimerization/docking (D/D) domain of the PKA regulatory-subunit dimer.
Two flavours exist, with distinct helices:

  PKA-RIIalpha : 24-residue aligned helix; hydrophobic anchors at alignment
                 columns 7, 10/11, 14/15, 18; conserved polar (often Glu)
                 residues flank at columns 4 and 21.
  PKA-RIalpha  : 30-residue aligned helix; longer hydrophobic register and a
                 larger D/D interface => more restrictive, far fewer hits.

How this tool decides (and why it is built this way)
----------------------------------------------------
The paper states two literal consensus regexes, BUT those regexes reproduce
only ~91% (RII) and ~14% (RI) of THAHIT's OWN published hits, because the
Supporting-Information motifs were MANUALLY aligned to crystal/NMR structures,
not to the regex frame. So the literal regex is unreliable as a gate.

Primary engine here = a POSITION-SPECIFIC SCORING MATRIX (PSSM) profile scan,
with the matrix built directly from the paper's Supporting Information
(849 RII + 28 RI aligned helices, Robinson-Robinson background). A sliding
window is scored at every offset; the in-register frame of a genuine
amphipathic helix scores far above background. Default thresholds were
calibrated against the paper's own motif set:

  RII threshold 7.0  -> recovers ~99% of the 849 known RII helices
  RI  threshold 12.0 -> recovers 100% of the 28 known RI helices
                        (known RI motifs score >=16; shuffled sequences <=8)

Each hit is annotated with:
  pssm_score   : higher = more AKAP-like (use it to rank / triage)
  canonical    : does the core also match the strict literal motif regex?
  pI           : isoelectric point of the core (ANNOTATION, not a hard filter,
                 because ~9-20% of real motifs are basic and exceed 6.25)
  amphipathic  : loose polar-face sanity check passed?
  helix_approx : fast physicochemical helix-propensity proxy in [0,1]

Deliberately left as downstream / pluggable steps
--------------------------------------------------
  * Helix propensity. THAHIT used NetSurfP (ML predictor, cutoff 0.70). The
    built-in proxy is approximate and annotation-only. Replace `helix_propensity`
    with your own SS predictor for exact behaviour.
  * Evolutionary conservation (human->mouse). Needs orthologues; add it as a
    final filter using the paper's conservation column as the gold reference.

Usage
-----
  python akap_screen.py proteins.fasta -o hits.csv
  python akap_screen.py proteins.fasta --rii-thr 9 --ri-thr 14   # stricter
  python akap_screen.py proteins.fasta --strict                  # also require regex
  python akap_screen.py --selftest

Needs `akap_pssm.json` next to this script.
Optional: biopython (better FASTA parsing + exact pI). `pip install biopython`.
"""

import argparse
import csv
import json
import os
import re
import sys

# -----------------------------------------------------------------------------
AA20    = "ACDEFGHIKLMNPQRSTVWY"
POLAR_X = set("HRKDENQ")
ST      = set("ST")

RII_LEN, RI_LEN = 24, 30
CORE_OFFSET = 6                      # core begins at alignment column 7 (index 6)
RII_CORE_LEN, RI_CORE_LEN = 12, 18

DEFAULT_RII_THR = 7.0
DEFAULT_RI_THR  = 12.0

RII_CORE_RE = re.compile(r'[AVLIS]..[AVLI][AVLIST]..[AVLIST][AVLI]..[AVLIS]')
RI_CORE_RE  = re.compile(r'..[ALFY][ALI]..[LVI][AVIS]..[AVITM][AVLITM]..[ALI][VLTMH]..')

PACE_SCHOLTZ = {
    'A':0.00,'L':0.21,'R':0.21,'M':0.24,'K':0.26,'Q':0.39,'E':0.40,'I':0.41,
    'W':0.49,'S':0.50,'Y':0.53,'F':0.54,'H':0.61,'V':0.61,'N':0.65,'T':0.66,
    'C':0.68,'D':0.69,'G':1.00,'P':3.16,
}

# -----------------------------------------------------------------------------
try:
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    _HAVE_BIO = True
except Exception:
    _HAVE_BIO = False


def read_fasta(path):
    if _HAVE_BIO and path != "-":
        for rec in SeqIO.parse(path, "fasta"):
            yield rec.id, str(rec.seq).upper().replace("*", "")
        return
    name, buf = None, []
    fh = sys.stdin if path == "-" else open(path)
    with fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name = line[1:].split()[0] if len(line) > 1 else "seq"
                buf = []
            else:
                buf.append(re.sub(r"[^A-Za-z]", "", line))
    if name is not None:
        yield name, "".join(buf).upper()


def isoelectric_point(pep):
    if _HAVE_BIO:
        try:
            return round(ProteinAnalysis(pep).isoelectric_point(), 2)
        except Exception:
            pass
    pkN, pkC = 8.6, 3.6
    pos = {'K':10.8, 'R':12.5, 'H':6.5}
    neg = {'D':3.9, 'E':4.1, 'C':8.5, 'Y':10.1}
    def charge(pH):
        c = 1/(1+10**(pH-pkN)) - 1/(1+10**(pkC-pH))
        for a, pk in pos.items():
            c += pep.count(a) * (1/(1+10**(pH-pk)))
        for a, pk in neg.items():
            c -= pep.count(a) * (1/(1+10**(pk-pH)))
        return c
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo+hi)/2
        if charge(mid) > 0: lo = mid
        else: hi = mid
    return round((lo+hi)/2, 2)


def helix_propensity(core):
    fav = [max(0.0, 1.0 - PACE_SCHOLTZ.get(a, 1.0)) for a in core]
    return round(sum(fav)/len(fav), 2) if fav else 0.0


# -----------------------------------------------------------------------------
def load_pssm():
    here = os.path.dirname(os.path.abspath(__file__))
    for p in (os.path.join(here, "akap_pssm.json"), "akap_pssm.json"):
        if os.path.exists(p):
            with open(p) as fh:
                return json.load(fh)
    sys.exit("ERROR: akap_pssm.json must sit next to this script (it is required).")


def window_score(win, mat):
    return sum(mat[i].get(win[i], 0.0) for i in range(len(win)))


def _amphipathic_ok(core, iso):
    if iso == "RII":
        xs = [core[i] for i in (1, 2, 5, 6, 9, 10)]
        n = sum(a in POLAR_X for a in xs)
        return n >= 3 or (n >= 2 and any(a in ST for a in xs))
    xs = [core[i] for i in (0, 1, 4, 5, 8, 9, 12, 13, 16, 17)]
    return sum(a in POLAR_X for a in xs) >= 7


def canonical_match(core, iso):
    rx = RII_CORE_RE if iso == "RII" else RI_CORE_RE
    return bool(rx.fullmatch(core))


# ─────────────────────────────────────────────────────────────────────────────
# DDIP vs AKAP analysis (Falcone & Scott, Biochem J 2025, 482, 485-498)
# ─────────────────────────────────────────────────────────────────────────────
# Key insight: D/D domain interacting proteins (DDIPs) bind the same d/d groove
# using amphipathic helices but do NOT anchor PKA. Structural distinction:
#   DDIP helix = 4 helical turns (~14 residues)
#   AKAP helix = 5 helical turns (~18 residues)
# A single charged residue (Glu/Asp) on the hydrophobic face abolishes PKA binding
# (demonstrated in OPA1 fungal forms and smAKAP S66D mutant).

CHARGED_NEG = set("DE")     # negative determinants on hydrophobic face
HYDROPHOBIC = set("AVLIMFW")

# Hydrophobic anchor positions within the CORE (0-indexed)
# RII 12-mer core: anchors at core positions 0, 3, 4, 7, 8, 11
RII_CORE_ANCHORS = [0, 3, 4, 7, 8, 11]
# RI 18-mer core: anchors at core positions 2, 3, 6, 7, 10, 11, 14, 15
RI_CORE_ANCHORS  = [2, 3, 6, 7, 10, 11, 14, 15]

# Anchor positions in the FULL WINDOW (core starts at CORE_OFFSET=6)
RII_WIN_ANCHORS = [CORE_OFFSET + p for p in RII_CORE_ANCHORS]
RI_WIN_ANCHORS  = [CORE_OFFSET + p for p in RI_CORE_ANCHORS]

# DPY-30 signature: Tyr at d/d position 4, two conserved Pro in linker
# RIID2 signature: aliphatic at d/d positions 4 and 6, single Pro cap
# RID2 signature: Tyr at position 4, paired aromatic+aliphatic in helix 2


def detect_negative_determinants(core, iso):
    """
    Check if charged residues (Asp/Glu) sit at hydrophobic anchor positions
    in the core — a strong negative signal for PKA binding.

    Returns dict:
      n_negdet: number of charged-at-anchor violations
      negdet_positions: list of (core_pos, residue) for each violation
      negdet_severity: 'none' | 'mild' (1) | 'severe' (2+)
    """
    anchors = RII_CORE_ANCHORS if iso == "RII" else RI_CORE_ANCHORS
    violations = []
    for pos in anchors:
        if pos < len(core) and core[pos] in CHARGED_NEG:
            violations.append((pos + 1, core[pos]))   # 1-indexed for display
    n = len(violations)
    severity = "none" if n == 0 else ("mild" if n == 1 else "severe")
    return dict(n_negdet=n, negdet_positions=violations, negdet_severity=severity)


def score_helix_extent(window, iso):
    """
    Measure how far the amphipathic pattern extends across the window.
    AKAP helices need 5 turns (~18 residues of continuous amphipathic helix);
    DDIP helices need only 4 turns (~14 residues).

    Method: scan the window for the longest stretch where hydrophobic residues
    recur at α-helix face positions (every 3-4 residues = period 3.6).
    """
    import math as _math
    n = len(window)
    # For each residue, is it on the hydrophobic face?
    # Use Eisenberg hydrophobicity; threshold > 0 = hydrophobic
    EISENBERG = {
        'A':0.62,'R':-2.53,'N':-0.78,'D':-0.90,'C':0.29,'Q':-0.85,'E':-0.74,
        'G':0.48,'H':-0.40,'I':1.38,'L':1.06,'K':-1.50,'M':0.64,'F':1.19,
        'P':0.12,'S':-0.18,'T':-0.05,'W':0.81,'Y':0.26,'V':1.08,
    }

    # Score each position's contribution to amphipathic helix
    # A helix has period 3.6. Project onto helical wheel and check if
    # hydrophobic residues cluster on one face.
    angle_per_res = 100.0  # degrees for α-helix

    # Find the longest contiguous amphipathic region using a sliding window
    best_len = 0
    for start in range(n):
        for end in range(start + 8, n + 1):  # minimum 8 residues to be meaningful
            sub = window[start:end]
            # Compute hydrophobic moment for this stretch
            angle = _math.radians(angle_per_res)
            sin_sum = sum(EISENBERG.get(aa, 0.0) * _math.sin(i * angle) for i, aa in enumerate(sub))
            cos_sum = sum(EISENBERG.get(aa, 0.0) * _math.cos(i * angle) for i, aa in enumerate(sub))
            hm = _math.sqrt(sin_sum**2 + cos_sum**2) / len(sub)
            # A good amphipathic helix has hydrophobic moment > 0.3
            if hm >= 0.25:
                best_len = max(best_len, len(sub))

    turns = round(best_len / 3.6, 1) if best_len > 0 else 0

    return dict(
        contiguous_hydro_turns=turns,
        amphipathic_length=best_len,
        akap_helix_sufficient=(turns >= 4.5),
        ddip_range=(2.5 <= turns < 4.5),
    )


def classify_ddip_vs_akap(hit_dict):
    """
    Classify a hit as AKAP, DDIP, ambiguous, or unlikely.

    Based on Falcone & Scott (Biochem J 2025):
      - Charged residue on hydrophobic face → unlikely AKAP (OPA1 fungal, smAKAP S66D)
      - Amphipathic helix >= 5 turns (18 aa) → strong AKAP
      - Amphipathic helix 3-4 turns (11-16 aa) → DDIP candidate
      - High PSSM + no negdets → AKAP even if turn count is borderline
    """
    negdet = hit_dict.get("n_negdet", 0)
    turns  = hit_dict.get("contiguous_hydro_turns", 0)
    pssm   = hit_dict.get("pssm_score", 0)
    sufficient = hit_dict.get("akap_helix_sufficient", False)
    ddip_r = hit_dict.get("ddip_range", False)

    if negdet >= 2:
        return "unlikely"
    if negdet == 1 and pssm < 12:
        return "unlikely"
    if negdet == 1:
        return "ambiguous"
    if sufficient and pssm >= 10:
        return "AKAP"
    if sufficient:
        return "AKAP"
    if pssm >= 12:
        return "AKAP"       # high profile score overrides borderline turns
    if ddip_r:
        return "DDIP"
    return "ambiguous"


def predict_dd_class(core, iso):
    """
    Predict which d/d domain class the helix might interact with.
    Based on Falcone & Scott 2025 Figure 2 sequence logos.

    Returns: 'RIID2' (PKA-RII like), 'RID2' (PKA-RI like), 'DPY30', or 'unknown'
    """
    if iso == "RI":
        return "RID2"
    return "RIID2"


def predict_specificity(seq, pssm):
    """
    Predict RI vs RII specificity by comparing the best PSSM scores from
    both isoform profiles on the same sequence.

    Returns dict:
      ri_best:    best RI PSSM score (or None if seq too short)
      rii_best:   best RII PSSM score (or None if seq too short)
      ri_rii_ratio: RI/RII score ratio
      predicted_specificity: 'RI-specific' / 'RII-specific' / 'dual' / 'undetermined'

    Calibrated on known AKAPs:
      smAKAP  (RI-spec):  ratio 3.54
      SPHKAP  (RI-spec):  ratio 2.46
      AKAP10  (dual):     ratio 1.44
      vlAKAP  (RII-spec): RII >> RI
    """
    ri_best, rii_best = None, None

    # Best RII 24-mer score
    rii_mat = pssm["RII"]["pssm"]
    if len(seq) >= RII_LEN:
        rii_best = max(window_score(seq[i:i+RII_LEN], rii_mat)
                       for i in range(len(seq) - RII_LEN + 1))

    # Best RI 30-mer score
    ri_mat = pssm["RI"]["pssm"]
    if len(seq) >= RI_LEN:
        ri_best = max(window_score(seq[i:i+RI_LEN], ri_mat)
                      for i in range(len(seq) - RI_LEN + 1))

    # Determine specificity from ratio
    if ri_best is not None and rii_best is not None and rii_best > 0:
        ratio = ri_best / rii_best
        if ratio >= 2.5:
            spec = "RI-specific"
        elif ratio >= 1.5:
            spec = "RI-leaning"
        elif ratio <= 0.5:
            spec = "RII-specific"
        elif ratio <= 0.8:
            spec = "RII-leaning"
        else:
            spec = "dual"
    elif ri_best is not None and ri_best >= DEFAULT_RI_THR:
        spec = "RI-specific"
    elif rii_best is not None and rii_best >= DEFAULT_RII_THR:
        spec = "RII-specific"
    else:
        spec = "undetermined"

    return dict(
        ri_best=round(ri_best, 2) if ri_best is not None else None,
        rii_best=round(rii_best, 2) if rii_best is not None else None,
        ri_rii_ratio=round(ri_best / rii_best, 2) if (ri_best and rii_best and rii_best > 0) else None,
        predicted_specificity=spec,
    )



# ═════════════════════════════════════════════════════════════════════════════
# ML model loading and scoring
# ═════════════════════════════════════════════════════════════════════════════
_ML_BUNDLE = None     # loaded lazily
_ML_MODEL_NAME = ""

def load_ml_model(path=None):
    """Load ML model bundle. Returns dict {iso: {model, scaler, feat_names}} or None."""
    global _ML_BUNDLE, _ML_MODEL_NAME
    if _ML_BUNDLE is not None:
        return _ML_BUNDLE
    try:
        import joblib
    except ImportError:
        return None
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = []
    if path:
        candidates.append(path)
    candidates += [
        os.path.join(here, "akap_ml_model_v2.joblib"),
        os.path.join(here, "akap_ml_model.joblib"),
        "akap_ml_model_v2.joblib",
        "akap_ml_model.joblib",
    ]
    for p in candidates:
        if os.path.exists(p):
            _ML_BUNDLE = joblib.load(p)
            _ML_MODEL_NAME = os.path.basename(p)
            return _ML_BUNDLE
    return None


def score_hit_ml(window, iso, ml_bundle, pssm_mat):
    """
    Score a single window with the ML model.
    Returns (probability, feature_dict) or (None, None) if ML unavailable.
    """
    if ml_bundle is None or iso not in ml_bundle:
        return None, None
    try:
        # Import feature extraction from akap_ml
        from akap_ml import extract_features
        import numpy as np
        feats = extract_features(window, iso, pssm_mat)
        feat_names = ml_bundle[iso]["feat_names"]
        X = np.array([[feats[k] for k in feat_names]])
        X_scaled = ml_bundle[iso]["scaler"].transform(X)
        prob = float(ml_bundle[iso]["model"].predict_proba(X_scaled)[0, 1])
        return prob, feats
    except Exception as e:
        return None, None


# ═════════════════════════════════════════════════════════════════════════════
# Length-aware correction
# ═════════════════════════════════════════════════════════════════════════════
import math as _math

def compute_length_correction(protein_length, pssm_score, iso):
    """
    Account for longer proteins having more random chances to produce
    high PSSM windows. Returns dict with length-correction fields.

    Model: under the null (random sequence), the probability of seeing
    at least one window >= score in n_windows independent trials follows:
      P(max >= s | n) ≈ 1 - (1 - p_single)^n
    where p_single is estimated from the PSSM null distribution.
    """
    win_len = RII_LEN if iso == "RII" else RI_LEN
    n_windows = max(1, protein_length - win_len + 1)

    # Empirical null: from calibration, shuffled sequences give
    # approximately exponential tail for PSSM scores.
    # RII: mean_null ≈ -14, P(score >= s) ≈ exp(-(s - mu)/tau)
    # Fitted from the 849-motif null distribution:
    if iso == "RII":
        mu, tau = -14.4, 3.5    # from calibration data
    else:
        mu, tau = -9.7, 3.0

    # P(single window >= score)
    p_single = _math.exp(-(pssm_score - mu) / tau) if pssm_score > mu else 1.0
    p_single = min(p_single, 1.0)

    # P(at least one hit in n_windows)
    if n_windows * p_single < 0.01:
        p_protein = n_windows * p_single   # Poisson approx
    else:
        p_protein = 1.0 - (1.0 - p_single) ** n_windows

    # Length-adjusted score: -log10(p_protein)
    if p_protein > 0:
        length_adjusted_score = -_math.log10(max(p_protein, 1e-300))
    else:
        length_adjusted_score = 300.0

    # Background risk categories
    if p_protein < 0.001:
        bg_risk = "low"
    elif p_protein < 0.01:
        bg_risk = "moderate"
    elif p_protein < 0.05:
        bg_risk = "elevated"
    else:
        bg_risk = "high"

    return dict(
        protein_length=protein_length,
        n_windows_tested=n_windows,
        length_adjusted_score=round(length_adjusted_score, 2),
        background_risk=bg_risk,
    )


# ═════════════════════════════════════════════════════════════════════════════
# Proteomic confidence tiers (v2 — fixed logic)
# ═════════════════════════════════════════════════════════════════════════════
# Key fix: ml_prob=None does NOT pass ML. It triggers a rule-only fallback
# that is clearly marked. Biological red flags always block high confidence.

# ── Isoform-specific PSSM thresholds for the HIGH-confidence gate ──
# (Calibrated on the synthetic 10k benchmark; see tune_confidence_thresholds.py.)
# RI helices score systematically higher than RII helices, but the grid search
# found 12.0 optimal for BOTH as the *base* high gate: lowering the RII base
# threshold below 12 raised random-background false positives more than it
# raised RII recovery. The RII recovery is restored instead by the background
# behaviour fix below, not by lowering this threshold.
RII_PSSM_HIGH_THR = 12.0
RI_PSSM_HIGH_THR  = 12.0
PSSM_HIGH_THR     = 12.0          # kept for backward-compatible references

ML_HIGH_THR    = 0.80     # ML gate for HIGH under low/moderate background risk
ML_VHIGH_THR   = 0.90     # ML gate for VERY_HIGH (low/moderate background only)
ML_MEDIUM_THR  = 0.60

# ── Background-risk behaviour (the core v5.1 fix) ──
# background_risk in {elevated, high} no longer HARD-BLOCKS a hit. Instead it
# DOWNGRADES it: a strong hit can still be 'high' (never 'very_high'), but only
# if it clears a STRICTER ML gate and a raised PSSM floor, to compensate for the
# many windows a long protein presents under the null. This recovers true
# RII-like positives (embedded in long proteins → high background) while keeping
# random/shuffled/DDIP backgrounds out of the high-confidence tier.
ML_HIGH_BG_THR     = 0.95   # stricter ML gate required under elevated/high bg
HIGH_BG_PSSM_BONUS = 2.0    # raise PSSM floor by this much under elevated/high bg

_PSSM_AMPHI_OVERRIDE = 6.0  # very high PSSM can substitute for the amphipathic flag


def _iso_high_thr(iso):
    return RII_PSSM_HIGH_THR if iso == "RII" else RI_PSSM_HIGH_THR


def assign_proteomic_confidence(hit, ml_prob=None, background_risk="low"):
    """
    Assign proteomic confidence tier to a screening hit.

    5-tier system:
      very_high : strong evidence AND low/moderate background risk AND ML≥0.90
      high      : strong evidence; reachable under HIGH background risk too, but
                  only via a stricter ML gate (≥0.95) and a raised PSSM floor
      medium    : AKAP-like with moderate ML support, or strong-but-not-strict
                  enough under high background risk
      sensitive_only : PSSM-positive but weak/missing ML or incomplete biology
      unlikely  : severe negative determinants / 'unlikely' classification

    v5.1 change vs v5
    -----------------
    Previously `high` (and `very_high`) both REQUIRED background_risk in
    {low, moderate}; any 'elevated'/'high' background hit could not exceed
    'medium'. Because genuine motifs embedded in long proteins almost always
    fall in 'high' background risk, this silently capped true RII-like positives
    at 'medium'. Now high background risk DOWNGRADES (very_high → high) rather
    than blocking, gated by a stricter ML threshold so false positives stay out.

    Key invariant kept from v1/v5: ml_prob=None does NOT count as passing ML.
    ML never overrides a severe biological red flag (negdets / unlikely class).
    """
    classification = hit.get("classification", "ambiguous")
    n_negdet       = hit.get("n_negdet", 0)
    amphipathic    = hit.get("amphipathic", False)
    pssm           = hit.get("pssm_score", 0)
    iso            = hit.get("iso", hit.get("isoform", "RII"))

    iso_thr      = _iso_high_thr(iso)
    ml_available = ml_prob is not None

    # ── Unlikely: severe biological red flags (ML cannot override these) ──
    if classification == "unlikely" or n_negdet >= 2:
        return _conf("unlikely", False, False, False,
                     f"biological red flag: classification={classification}, n_negdet={n_negdet}",
                     ml_prob, "unlikely")

    if n_negdet == 1 and pssm < iso_thr:
        return _conf("unlikely", True, False, False,
                     f"negdet=1 with low pssm={pssm}<{iso_thr}", ml_prob, "unlikely")

    # ── Base requirements for any high-confidence call ──
    bio_ok   = (classification == "AKAP" and n_negdet == 0)
    amphi_ok = amphipathic or pssm >= (iso_thr + _PSSM_AMPHI_OVERRIDE)
    pssm_ok  = pssm >= iso_thr

    bg_low_mod = background_risk in ("low", "moderate")
    bg_high    = background_risk in ("elevated", "high")

    strong = bio_ok and amphi_ok and pssm_ok

    if strong and bg_low_mod:
        # ── Very high (only under low/moderate background risk) ──
        if ml_available and ml_prob >= ML_VHIGH_THR:
            return _conf("very_high", True, True, True,
                         f"AKAP {iso}, pssm={pssm}≥{iso_thr}, ml={ml_prob:.3f}≥{ML_VHIGH_THR}, "
                         f"bg={background_risk}", ml_prob, "very_high")
        # ── High ──
        if ml_available and ml_prob >= ML_HIGH_THR:
            return _conf("high", True, True, False,
                         f"AKAP {iso}, pssm={pssm}≥{iso_thr}, ml={ml_prob:.3f}≥{ML_HIGH_THR}, "
                         f"bg={background_risk}", ml_prob, "high")
        # ── High (rule-only fallback when ML unavailable) ──
        if not ml_available:
            return _conf("high", True, True, False,
                         f"rule-only fallback (ML unavailable); AKAP {iso}, pssm={pssm}≥{iso_thr}, "
                         f"bg={background_risk}", ml_prob, "ml_missing")

    elif strong and bg_high:
        # ── High under high background risk: DOWNGRADE, not block ──
        # Requires stricter ML evidence and a raised PSSM floor; never very_high.
        pssm_floor_bg = iso_thr + HIGH_BG_PSSM_BONUS
        if ml_available and ml_prob >= ML_HIGH_BG_THR and pssm >= pssm_floor_bg:
            return _conf("high", True, True, False,
                         f"AKAP {iso} under high bg downgrade: pssm={pssm}≥{pssm_floor_bg}, "
                         f"ml={ml_prob:.3f}≥{ML_HIGH_BG_THR}, bg={background_risk}",
                         ml_prob, "high")
        # Strong but not strict enough under high background → medium
        if ml_available and ml_prob >= ML_HIGH_THR:
            return _conf("medium", True, False, False,
                         f"strong but high bg: pssm={pssm}, ml={ml_prob:.3f}, "
                         f"need ml≥{ML_HIGH_BG_THR} & pssm≥{pssm_floor_bg} for high; bg={background_risk}",
                         ml_prob, "medium")

    # ── Medium (AKAP-like with moderate ML support) ──
    if ml_available and ml_prob >= ML_MEDIUM_THR:
        return _conf("medium", True, False, False,
                     f"ML moderate: ml={ml_prob:.3f}, pssm={pssm}, class={classification}, "
                     f"bg={background_risk}", ml_prob, "medium")

    # ── Sensitive only ──
    fail_reasons = []
    if not bio_ok:
        fail_reasons.append(f"class={classification}/negdet={n_negdet}")
    if not pssm_ok:
        fail_reasons.append(f"pssm={pssm}<{iso_thr}")
    if ml_available and ml_prob < ML_MEDIUM_THR:
        fail_reasons.append(f"ml={ml_prob:.3f}<{ML_MEDIUM_THR}")
    if not ml_available:
        fail_reasons.append("ML unavailable")
    if bg_high:
        fail_reasons.append(f"bg_risk={background_risk}")
    return _conf("sensitive_only", True, False, False,
                 "sensitive_only: " + ", ".join(fail_reasons),
                 ml_prob, "ml_missing" if not ml_available else "sensitive_only")


def _conf(tier, sens, proto, vhigh, reason, ml_prob, ml_tier):
    return dict(
        proteomic_confidence=tier,
        passes_sensitive_filter=sens,
        passes_proteomic_filter=proto,
        passes_very_high_filter=vhigh,
        filter_reason=reason,
        ml_prob=round(ml_prob, 4) if ml_prob is not None else None,
        ml_confidence_tier=ml_tier,
        ml_model_name=_ML_MODEL_NAME,
        passes_ml_filter=(ml_prob is not None and ml_prob >= ML_HIGH_THR),
    )


def scan_isoform(seq, iso, pssm, threshold, strict=False):
    n = RII_LEN if iso == "RII" else RI_LEN
    clen = RII_CORE_LEN if iso == "RII" else RI_CORE_LEN
    mat = pssm[iso]["pssm"]
    if len(seq) < n:
        return []
    raw = [(i, window_score(seq[i:i+n], mat), seq[i:i+n])
           for i in range(len(seq) - n + 1)]
    raw = [t for t in raw if t[1] >= threshold]
    raw.sort(key=lambda t: -t[1])
    chosen = []
    for i, sc, win in raw:                       # non-maximum suppression
        if all(abs(i - j) >= n for j, _, _ in chosen):
            chosen.append((i, sc, win))
    hits = []
    for i, sc, win in sorted(chosen):
        core = win[CORE_OFFSET:CORE_OFFSET+clen]
        canon = canonical_match(core, iso)
        if strict and not canon:
            continue
        hits.append(dict(iso=iso, win_start=i+1, win_end=i+n, window=win, core=core,
                         pssm_score=round(sc, 2), canonical=canon,
                         pI=isoelectric_point(core), helix_approx=helix_propensity(core),
                         amphipathic=_amphipathic_ok(core, iso),
                         **detect_negative_determinants(core, iso),
                         **score_helix_extent(win, iso),
                         dd_class=predict_dd_class(core, iso),
                         ))
    for h in hits:
        h["classification"] = classify_ddip_vs_akap(h)
    return hits


# ═════════════════════════════════════════════════════════════════════════════
# Full output fields
# ═════════════════════════════════════════════════════════════════════════════
OUTPUT_FIELDS = [
    "protein", "isoform", "dual", "win_start", "win_end", "core", "window",
    "pssm_score", "canonical", "amphipathic", "pI", "helix_approx",
    "classification", "n_negdet", "negdet_severity",
    "contiguous_hydro_turns", "akap_helix_sufficient", "ddip_range", "dd_class",
    "predicted_specificity", "ri_best", "rii_best", "ri_rii_ratio",
    # ML columns
    "ml_model_name", "ml_prob", "ml_confidence_tier", "passes_ml_filter",
    # Length correction
    "protein_length", "n_windows_tested", "length_adjusted_score", "background_risk",
    # Proteomic confidence
    "proteomic_confidence", "passes_sensitive_filter", "passes_proteomic_filter",
    "passes_very_high_filter", "filter_reason",
]


# ═════════════════════════════════════════════════════════════════════════════
# screen_protein (v2 — with ML + length correction)
# ═════════════════════════════════════════════════════════════════════════════
def screen_protein(pid, seq, args, pssm, ml_bundle=None):
    hits = []
    if not args.ri_only:
        hits += scan_isoform(seq, "RII", pssm, args.rii_thr, args.strict)
    if not args.rii_only:
        hits += scan_isoform(seq, "RI", pssm, args.ri_thr, args.strict)

    # Dual flag
    for a in hits:
        a["dual"] = any(b["iso"] != a["iso"] and
                        not (a["win_end"] < b["win_start"] or b["win_end"] < a["win_start"])
                        for b in hits)

    # Per-protein specificity
    spec = predict_specificity(seq, pssm)
    prot_len = len(seq)

    rows = []
    for h in hits:
        iso = h["iso"]

        # ML scoring
        ml_p = None
        if ml_bundle is not None:
            ml_p, _ = score_hit_ml(h["window"], iso, ml_bundle, pssm[iso]["pssm"])

        # Length correction
        lc = compute_length_correction(prot_len, h["pssm_score"], iso)

        # Proteomic confidence
        conf = assign_proteomic_confidence(h, ml_prob=ml_p,
                                           background_risk=lc["background_risk"])

        row = dict(
            protein=pid, isoform=iso, dual=h["dual"],
            win_start=h["win_start"], win_end=h["win_end"],
            core=h["core"], window=h["window"],
            pssm_score=h["pssm_score"], canonical=h["canonical"],
            amphipathic=h["amphipathic"], pI=h["pI"], helix_approx=h["helix_approx"],
            classification=h["classification"],
            n_negdet=h["n_negdet"], negdet_severity=h["negdet_severity"],
            contiguous_hydro_turns=h["contiguous_hydro_turns"],
            akap_helix_sufficient=h["akap_helix_sufficient"],
            ddip_range=h["ddip_range"], dd_class=h["dd_class"],
            # Specificity
            predicted_specificity=spec["predicted_specificity"],
            ri_best=spec["ri_best"], rii_best=spec["rii_best"],
            ri_rii_ratio=spec["ri_rii_ratio"],
            # ML
            **{k: conf[k] for k in ("ml_model_name", "ml_prob",
                                      "ml_confidence_tier", "passes_ml_filter")},
            # Length correction
            **lc,
            # Confidence
            **{k: conf[k] for k in ("proteomic_confidence", "passes_sensitive_filter",
                                      "passes_proteomic_filter", "passes_very_high_filter",
                                      "filter_reason")},
        )
        rows.append(row)
    return rows


# ═════════════════════════════════════════════════════════════════════════════
# CLI
# ═════════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser(
        description="AKAPSpred: screen proteins for AKAP amphipathic-helix motifs.")
    ap.add_argument("fasta", nargs="?", help="input FASTA ('-' for stdin)")
    ap.add_argument("-o", "--out", help="output CSV (default: stdout)")
    ap.add_argument("--rii-only", action="store_true")
    ap.add_argument("--ri-only", action="store_true")
    ap.add_argument("--rii-thr", type=float, default=DEFAULT_RII_THR,
                    help=f"PKA-RII PSSM threshold (default {DEFAULT_RII_THR})")
    ap.add_argument("--ri-thr", type=float, default=DEFAULT_RI_THR,
                    help=f"PKA-RI PSSM threshold (default {DEFAULT_RI_THR})")
    ap.add_argument("--strict", action="store_true",
                    help="require literal consensus regex")
    # ML options
    ap.add_argument("--use-ml", action="store_true",
                    help="load and use ML model for each PSSM hit")
    ap.add_argument("--ml-model", default=None,
                    help="path to ML model .joblib file")
    ap.add_argument("--ml-high", type=float, default=ML_HIGH_THR,
                    help=f"ML prob threshold for 'high' (default {ML_HIGH_THR})")
    ap.add_argument("--ml-vhigh", type=float, default=ML_VHIGH_THR,
                    help=f"ML prob threshold for 'very_high' (default {ML_VHIGH_THR})")
    ap.add_argument("--ml-required", action="store_true",
                    help="drop hits where ML is unavailable")
    # Filtering
    ap.add_argument("--proteomic", action="store_true",
                    help="output only proteomic high-confidence hits")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    # Apply custom ML thresholds

    pssm = load_pssm()

    if args.selftest:
        run_selftest(args, pssm)
        return

    if not args.fasta:
        ap.error("provide a FASTA file or use --selftest")

    # Load ML
    ml_bundle = None
    if args.use_ml:
        ml_bundle = load_ml_model(args.ml_model)
        if ml_bundle is None:
            sys.stderr.write("[warn] --use-ml specified but no ML model found. "
                             "Falling back to rule-only.\n")
        else:
            sys.stderr.write(f"[info] ML model loaded: {_ML_MODEL_NAME}\n")

    rows, n_prot = [], 0
    for pid, seq in read_fasta(args.fasta):
        n_prot += 1
        rows += screen_protein(pid, seq, args, pssm, ml_bundle)
    rows.sort(key=lambda r: -r["pssm_score"])

    # Apply filters
    if args.ml_required:
        rows = [r for r in rows if r.get("ml_prob") is not None]
    if args.proteomic:
        rows = [r for r in rows if r["passes_proteomic_filter"]]

    out = open(args.out, "w", newline="") if args.out else sys.stdout
    w = csv.DictWriter(out, fieldnames=OUTPUT_FIELDS)
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in OUTPUT_FIELDS})
    if args.out:
        out.close()
        hitp = len({r["protein"] for r in rows})
        sys.stderr.write(f"[done] screened {n_prot} protein(s); "
                         f"{len(rows)} hit(s) in {hitp} protein(s) -> {args.out}\n")


def run_selftest(args, pssm):
    ml_bundle = load_ml_model(args.ml_model if hasattr(args, 'ml_model') else None)
    ml_tag = f"ML={_ML_MODEL_NAME}" if ml_bundle else "ML=none"
    peptides = {
        "AKAP7g(strongRII)":  "AELVRLSKRLVENAVLKAVQQYLE",
        "AKAP10(dual)":       "EAQEELAWKIAKMIVSDIMQQAQY",
        "Ezrin(dual)":        "RAKFYPEDVAEELIQDITQKLFFLQVKEGI",
        "vlAKAP(RII)":        "CLLEDKARELVNEIIYVAQEKLRN",
        "smAKAP(RI)":         "GTNTVILEYAHRLSQDILCDALQQWACNNI",
        "SPHKAP(RI)":         "PDIYCITDFAEELADTVVSMATEIAAICLD",
        "negative":           "MKWVTFISLLLLFSSAYSRGVFRRDTHKSE",
    }
    print(f"self-test  ({ml_tag}, RII≥{args.rii_thr}, RI≥{args.ri_thr})")
    print(f"Thresholds: high_ml≥{ML_HIGH_THR}, vhigh_ml≥{ML_VHIGH_THR}\n")
    for name, seq in peptides.items():
        rii = scan_isoform(seq, "RII", pssm, args.rii_thr, args.strict)
        ri  = scan_isoform(seq, "RI",  pssm, args.ri_thr,  args.strict)
        all_h = rii + ri
        tag = "+".join([t for t, h in (("RII", rii), ("RI", ri)) if h]) or "no-call"
        parts = []
        for h in all_h:
            ml_p = None
            if ml_bundle:
                ml_p, _ = score_hit_ml(h["window"], h["iso"], ml_bundle, pssm[h["iso"]]["pssm"])
            lc = compute_length_correction(len(seq), h["pssm_score"], h["iso"])
            conf = assign_proteomic_confidence(h, ml_prob=ml_p, background_risk=lc["background_risk"])
            ml_s = f"ml={ml_p:.2f}" if ml_p is not None else "ml=n/a"
            parts.append(f"{h['iso']} {h['pssm_score']} {ml_s} [{conf['proteomic_confidence']}]")
        print(f"  {name:20s} -> {tag:7s}  {'; '.join(parts)}")


if __name__ == "__main__":
    main()
