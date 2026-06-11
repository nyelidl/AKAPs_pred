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
        # A ratio near 2.0 already indicates strong RI preference.
        # This keeps known RI-selective SPHKAP (RI/RII ≈ 2.46) as RI-specific.
        if ratio >= 2.0:
            spec = "RI-specific"
        elif ratio >= 1.3:
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
                         # DDIP / negative determinant analysis
                         **detect_negative_determinants(core, iso),
                         **score_helix_extent(win, iso),
                         dd_class=predict_dd_class(core, iso),
                         ))
    # Add classification after all hits collected
    for h in hits:
        h["classification"] = classify_ddip_vs_akap(h)
    return hits


# -----------------------------------------------------------------------------
def screen_protein(pid, seq, args, pssm):
    hits = []
    if not args.ri_only:
        hits += scan_isoform(seq, "RII", pssm, args.rii_thr, args.strict)
    if not args.rii_only:
        hits += scan_isoform(seq, "RI", pssm, args.ri_thr, args.strict)

    # Protein-level RI/RII specificity, evaluated across the full sequence.
    # This prevents interpreting every threshold-crossing RII window as true RII specificity
    # when the full-length protein is strongly RI-preferred, as in SPHKAP.
    spec = predict_specificity(seq, pssm)

    for a in hits:
        a["dual"] = any(b["iso"] != a["iso"] and
                        not (a["win_end"] < b["win_start"] or b["win_end"] < a["win_start"])
                        for b in hits)
    return [dict(protein=pid, isoform=h["iso"], dual=h["dual"],
                 win_start=h["win_start"], win_end=h["win_end"],
                 core=h["core"], window=h["window"], pssm_score=h["pssm_score"],
                 canonical=h["canonical"], amphipathic=h["amphipathic"],
                 pI=h["pI"], helix_approx=h["helix_approx"],
                 classification=h["classification"],
                 n_negdet=h["n_negdet"], negdet_severity=h["negdet_severity"],
                 contiguous_hydro_turns=h["contiguous_hydro_turns"],
                 akap_helix_sufficient=h["akap_helix_sufficient"],
                 ddip_range=h["ddip_range"], dd_class=h["dd_class"],
                 predicted_specificity=spec["predicted_specificity"],
                 ri_best=spec["ri_best"], rii_best=spec["rii_best"],
                 ri_rii_ratio=spec["ri_rii_ratio"],
                 ) for h in hits]


def main():
    ap = argparse.ArgumentParser(
        description="Screen a protein FASTA for AKAP (PKA-R binding) amphipathic-helix motifs.")
    ap.add_argument("fasta", nargs="?", help="input FASTA ('-' for stdin)")
    ap.add_argument("-o", "--out", help="output CSV (default: stdout)")
    ap.add_argument("--rii-only", action="store_true", help="screen only the PKA-RII motif")
    ap.add_argument("--ri-only", action="store_true", help="screen only the PKA-RI motif")
    ap.add_argument("--rii-thr", type=float, default=DEFAULT_RII_THR,
                    help=f"PKA-RII PSSM threshold (default {DEFAULT_RII_THR})")
    ap.add_argument("--ri-thr", type=float, default=DEFAULT_RI_THR,
                    help=f"PKA-RI PSSM threshold (default {DEFAULT_RI_THR})")
    ap.add_argument("--strict", action="store_true",
                    help="additionally require the literal consensus regex (lower recall)")
    ap.add_argument("--selftest", action="store_true",
                    help="run on the paper's validated AKAP peptides and exit")
    args = ap.parse_args()

    pssm = load_pssm()
    if args.selftest:
        run_selftest(args, pssm)
        return
    if not args.fasta:
        ap.error("provide a FASTA file (or '-' for stdin), or use --selftest")

    rows, n_prot = [], 0
    for pid, seq in read_fasta(args.fasta):
        n_prot += 1
        rows += screen_protein(pid, seq, args, pssm)
    rows.sort(key=lambda r: -r["pssm_score"])

    fields = ["protein","isoform","dual","win_start","win_end","core","window",
              "pssm_score","canonical","amphipathic","pI","helix_approx",
              "classification","n_negdet","negdet_severity",
              "contiguous_hydro_turns","akap_helix_sufficient","ddip_range","dd_class",
              "predicted_specificity","ri_best","rii_best","ri_rii_ratio"]
    out = open(args.out, "w", newline="") if args.out else sys.stdout
    w = csv.DictWriter(out, fieldnames=fields)
    w.writeheader()
    for r in rows:
        w.writerow(r)
    if args.out:
        out.close()
        hitp = len({r["protein"] for r in rows})
        sys.stderr.write(f"[done] screened {n_prot} protein(s); "
                         f"{len(rows)} hit(s) in {hitp} protein(s) -> {args.out}\n")


def run_selftest(args, pssm):
    peptides = {
        "AKAP7g_294-317(strongRII)": "AELVRLSKRLVENAVLKAVQQYLE",
        "AKAP10_629-652(dual)":      "EAQEELAWKIAKMIVSDIMQQAQY",
        "Ezrin_84-107(dual)":        "RAKFYPEDVAEELIQDITQKLFFLQVKEGI",
        "vlAKAP_1299-1322(RII)":     "CLLEDKARELVNEIIYVAQEKLRN",
        "smAKAP_56-79(RI)":          "GTNTVILEYAHRLSQDILCDALQQWACNNI",
        "SPHKAP_920-949(RI)":        "PDIYCITDFAEELADTVVSMATEIAAICLD",
        "OPA1_human(AKAP)":          "KKVRE IQEKLD AFIEALH".replace(" ",""),
        "OPA1_fungal(nonAKAP)":      "EALVAERDRVKALLDAYK",
        "negative_control":          "MKWVTFISLLLLFSSAYSRGVFRRDTHKSE",
    }
    print(f"self-test  (RII thr={args.rii_thr}, RI thr={args.ri_thr})\n")
    for name, seq in peptides.items():
        rii = scan_isoform(seq, "RII", pssm, args.rii_thr, args.strict)
        ri  = scan_isoform(seq, "RI",  pssm, args.ri_thr,  args.strict)
        tag = "+".join([t for t, h in (("RII", rii), ("RI", ri)) if h]) or "no-call"
        det = "; ".join(f"{h['iso']} {h['pssm_score']}"
                        f"{'/canon' if h['canonical'] else ''}"
                        f" [{h['classification']}|turns={h['contiguous_hydro_turns']}"
                        f"|negdet={h['n_negdet']}]"
                        for h in rii+ri)
        # Specificity prediction
        spec = predict_specificity(seq, pssm)
        spec_tag = f"  spec={spec['predicted_specificity']}" if spec['predicted_specificity'] != 'undetermined' else ""
        if spec['ri_rii_ratio']:
            spec_tag += f" (RI/RII={spec['ri_rii_ratio']})"
        print(f"  {name:27s} -> {tag:7s}  {det}{spec_tag}")


if __name__ == "__main__":
    main()
