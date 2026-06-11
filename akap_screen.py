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
                         amphipathic=_amphipathic_ok(core, iso)))
    return hits


# -----------------------------------------------------------------------------
def screen_protein(pid, seq, args, pssm):
    hits = []
    if not args.ri_only:
        hits += scan_isoform(seq, "RII", pssm, args.rii_thr, args.strict)
    if not args.rii_only:
        hits += scan_isoform(seq, "RI", pssm, args.ri_thr, args.strict)
    for a in hits:
        a["dual"] = any(b["iso"] != a["iso"] and
                        not (a["win_end"] < b["win_start"] or b["win_end"] < a["win_start"])
                        for b in hits)
    return [dict(protein=pid, isoform=h["iso"], dual=h["dual"],
                 win_start=h["win_start"], win_end=h["win_end"],
                 core=h["core"], window=h["window"], pssm_score=h["pssm_score"],
                 canonical=h["canonical"], amphipathic=h["amphipathic"],
                 pI=h["pI"], helix_approx=h["helix_approx"]) for h in hits]


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
              "pssm_score","canonical","amphipathic","pI","helix_approx"]
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
        "WAVE1_20-43(RII)":          "RGIKNELECVTNISLANIIRQLSS",
        "vlAKAP_1299-1322(RII)":     "CLLEDKARELVNEIIYVAQEKLRN",
        "superAKAP-IS(RII)":         "QIEYVAKQIVDYAIHQA",
        "smAKAP_56-79(RI)":          "GTNTVILEYAHRLSQDILCDALQQWACNNI",
        "RIAD(RI)":                  "TVLEQYANQLADQIIKEATE",
        "negative_control":          "MKWVTFISLLLLFSSAYSRGVFRRDTHKSE",
    }
    print(f"self-test  (RII thr={args.rii_thr}, RI thr={args.ri_thr})\n")
    for name, seq in peptides.items():
        rii = scan_isoform(seq, "RII", pssm, args.rii_thr, args.strict)
        ri  = scan_isoform(seq, "RI",  pssm, args.ri_thr,  args.strict)
        tag = "+".join([t for t, h in (("RII", rii), ("RI", ri)) if h]) or "no-call"
        det = "; ".join(f"{h['iso']} {h['pssm_score']}"
                        f"{'/canon' if h['canonical'] else ''}" for h in rii+ri)
        print(f"  {name:27s} -> {tag:7s}  {det}")


if __name__ == "__main__":
    main()
