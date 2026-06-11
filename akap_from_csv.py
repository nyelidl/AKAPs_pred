#!/usr/bin/env python3
"""
akap_from_csv.py
================
Wrapper that reads a CSV of proteins, obtains sequences (from the CSV itself
or by fetching from UniProt), then runs the AKAP motif screen.

Input CSV format (flexible — column names are auto-detected):
-------------------------------------------------------------
The script looks for columns whose header contains these keywords (case-insensitive):

  sequence / seq        -> amino-acid sequence (used directly if present)
  uniprot / accession   -> UniProt accession ID (used to fetch if no sequence)
  protein / name / gene -> display name (optional; falls back to UniProt ID)

Any column order works. Extra columns are ignored.

Examples of valid CSV headers:
  protein_name,uniprot_id,sequence
  Name,Accession
  gene,seq
  UniProt                              (single column is fine)

Workflow
--------
1. For rows that already have a sequence -> use it.
2. For rows without a sequence but with a UniProt ID -> fetch from UniProt API.
3. Write a combined FASTA file.
4. Run akap_screen.py and produce the final hits CSV.

Usage
-----
  python akap_from_csv.py my_proteins.csv -o hits.csv
  python akap_from_csv.py my_proteins.csv -o hits.csv --rii-thr 9
  python akap_from_csv.py my_proteins.csv --keep-fasta         # also keep the FASTA

Requirements
------------
  - akap_screen.py + akap_pssm.json (same directory or specify --script-dir)
  - requests (pip install requests) — only needed if fetching from UniProt
  - pandas   (pip install pandas)
"""

import argparse
import csv
import os
import re
import sys
import time
import tempfile

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas is required. Install with: pip install pandas")


# ---------------------------------------------------------------------------
# Column detection
# ---------------------------------------------------------------------------
def _find_col(columns, keywords):
    """Return the first column name that contains any of the keywords (case-insensitive)."""
    for col in columns:
        low = col.lower().strip()
        for kw in keywords:
            if kw in low:
                return col
    return None


def detect_columns(df):
    """Auto-detect which columns hold name, uniprot ID, and sequence."""
    cols = list(df.columns)
    seq_col = _find_col(cols, ["sequence", "seq"])
    uni_col = _find_col(cols, ["uniprot", "accession", "acc_id", "entry"])
    name_col = _find_col(cols, ["protein", "name", "gene", "description", "label"])

    # If only one column and none matched, guess by content
    if len(cols) == 1 and not any([seq_col, uni_col, name_col]):
        sample = str(df.iloc[0, 0]).strip()
        if re.match(r'^[A-Za-z]{10,}$', sample):
            seq_col = cols[0]
        elif re.match(r'^[A-Z][0-9A-Z]{5}$', sample) or re.match(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$', sample):
            uni_col = cols[0]
        else:
            name_col = cols[0]

    return seq_col, uni_col, name_col


# ---------------------------------------------------------------------------
# UniProt fetcher
# ---------------------------------------------------------------------------
def fetch_uniprot_sequence(accession, retries=2, timeout=10):
    """Fetch a single protein sequence from UniProt REST API. Returns (sequence, error)."""
    try:
        import requests
    except ImportError:
        return None, "requests library not installed (pip install requests)"

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(retries + 1):
        try:
            resp = requests.get(url, timeout=timeout)
            if resp.status_code == 200:
                lines = resp.text.strip().split("\n")
                seq = "".join(l.strip() for l in lines if not l.startswith(">"))
                if seq:
                    return seq.upper(), None
                return None, "empty sequence returned"
            elif resp.status_code == 404:
                return None, f"not found (404)"
            else:
                if attempt < retries:
                    time.sleep(1)
                    continue
                return None, f"HTTP {resp.status_code}"
        except Exception as e:
            if attempt < retries:
                time.sleep(1)
                continue
            return None, str(e)
    return None, "max retries exceeded"


def fetch_uniprot_batch(accessions, verbose=True):
    """Fetch sequences for a list of UniProt accessions. Returns dict {accession: sequence} and list of errors."""
    results = {}
    errors = []
    total = len(accessions)
    for i, acc in enumerate(accessions, 1):
        acc = acc.strip()
        if not acc:
            continue
        if verbose:
            print(f"  [{i}/{total}] Fetching {acc}...", end=" ", flush=True)
        seq, err = fetch_uniprot_sequence(acc)
        if seq:
            results[acc] = seq
            if verbose:
                print(f"OK ({len(seq)} aa)")
        else:
            errors.append((acc, err))
            if verbose:
                print(f"FAILED: {err}")
        # Be polite to the API
        if i < total:
            time.sleep(0.2)
    return results, errors


# ---------------------------------------------------------------------------
# FASTA writer
# ---------------------------------------------------------------------------
def write_fasta(entries, path):
    """Write [(id, sequence), ...] to a FASTA file."""
    with open(path, "w") as f:
        for pid, seq in entries:
            # Clean the sequence
            seq_clean = re.sub(r'[^A-Za-z]', '', seq).upper()
            f.write(f">{pid}\n")
            # Wrap at 80 chars
            for i in range(0, len(seq_clean), 80):
                f.write(seq_clean[i:i+80] + "\n")
    return len(entries)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Read a CSV of proteins, fetch sequences if needed, and screen for AKAP motifs.")
    ap.add_argument("csv_file", help="input CSV file (protein names / UniProt IDs / sequences)")
    ap.add_argument("-o", "--out", default="akap_hits.csv",
                    help="output CSV with AKAP hits (default: akap_hits.csv)")
    ap.add_argument("--keep-fasta", action="store_true",
                    help="keep the intermediate FASTA file (saved as <csv_name>.fasta)")
    ap.add_argument("--fasta-out", default=None,
                    help="path for the intermediate FASTA (default: auto)")
    ap.add_argument("--rii-thr", type=float, default=7.0,
                    help="PKA-RII PSSM threshold (default 7.0)")
    ap.add_argument("--ri-thr", type=float, default=12.0,
                    help="PKA-RI PSSM threshold (default 12.0)")
    ap.add_argument("--rii-only", action="store_true", help="screen only PKA-RII")
    ap.add_argument("--ri-only", action="store_true", help="screen only PKA-RI")
    ap.add_argument("--strict", action="store_true",
                    help="also require the literal consensus regex")
    ap.add_argument("--script-dir", default=None,
                    help="directory containing akap_screen.py and akap_pssm.json")
    args = ap.parse_args()

    # --- Locate akap_screen ---
    script_dir = args.script_dir or os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, script_dir)
    try:
        import akap_screen as A
    except ImportError:
        sys.exit(f"ERROR: cannot import akap_screen.py from {script_dir}\n"
                 f"       Place akap_screen.py + akap_pssm.json there, or use --script-dir.")
    pssm = A.load_pssm()

    # --- Read CSV ---
    print(f"\nReading {args.csv_file}...")
    # Try to detect separator
    with open(args.csv_file) as f:
        sample = f.read(4096)
    sep = ","
    if sample.count("\t") > sample.count(","):
        sep = "\t"
    elif sample.count(";") > sample.count(","):
        sep = ";"

    df = pd.read_csv(args.csv_file, sep=sep, dtype=str).fillna("")
    # Strip whitespace from headers and values
    df.columns = [c.strip() for c in df.columns]
    df = df.apply(lambda col: col.str.strip())

    print(f"  {len(df)} rows, columns: {list(df.columns)}")

    seq_col, uni_col, name_col = detect_columns(df)
    print(f"  Detected: name='{name_col}', uniprot='{uni_col}', sequence='{seq_col}'")

    # --- Build protein entries ---
    entries = []            # [(display_id, sequence), ...]
    need_fetch = []         # [(row_idx, accession), ...]

    for idx, row in df.iterrows():
        # Get display name
        name = ""
        if name_col:
            name = row[name_col]
        if not name and uni_col:
            name = row[uni_col]
        if not name:
            name = f"protein_{idx+1}"
        # Clean name for FASTA header (no spaces/special chars)
        display = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)

        # Check for sequence
        seq = ""
        if seq_col:
            raw = row[seq_col]
            seq = re.sub(r'[^A-Za-z]', '', raw).upper()

        if len(seq) >= 20:
            entries.append((display, seq))
        elif uni_col and row[uni_col]:
            need_fetch.append((len(entries), display, row[uni_col]))
            entries.append((display, ""))   # placeholder
        else:
            print(f"  [warn] Row {idx+1} ({name}): no sequence and no UniProt ID — skipped.")

    # --- Fetch from UniProt if needed ---
    if need_fetch:
        print(f"\n  {len(need_fetch)} protein(s) need sequence from UniProt API...")
        accessions = [acc for _, _, acc in need_fetch]
        fetched, errors = fetch_uniprot_batch(accessions)

        for eidx, display, acc in need_fetch:
            if acc in fetched:
                entries[eidx] = (display, fetched[acc])
            else:
                print(f"  [warn] Could not fetch {acc} ({display}) — will be skipped.")

    # Filter out empty sequences
    entries = [(pid, seq) for pid, seq in entries if len(seq) >= 20]
    print(f"\n  {len(entries)} protein(s) with valid sequences ready for screening.")

    if not entries:
        sys.exit("ERROR: no valid protein sequences to screen.")

    # --- Write FASTA ---
    if args.fasta_out:
        fasta_path = args.fasta_out
    elif args.keep_fasta:
        base = os.path.splitext(os.path.basename(args.csv_file))[0]
        fasta_path = base + ".fasta"
    else:
        tmp = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w")
        fasta_path = tmp.name
        tmp.close()

    n = write_fasta(entries, fasta_path)
    print(f"  Wrote {n} sequences to {fasta_path}")

    # --- Run the screen ---
    print(f"\nRunning AKAP screen (RII thr={args.rii_thr}, RI thr={args.ri_thr})...\n")

    all_rows = []
    for pid, seq in entries:
        # Build a minimal args-like object for screen_protein
        hits = []
        if not args.ri_only:
            hits += A.scan_isoform(seq, "RII", pssm, args.rii_thr, args.strict)
        if not args.rii_only:
            hits += A.scan_isoform(seq, "RI", pssm, args.ri_thr, args.strict)
        # dual flag
        for a in hits:
            a["dual"] = any(b["iso"] != a["iso"] and
                            not (a["win_end"] < b["win_start"] or b["win_end"] < a["win_start"])
                            for b in hits)
        for h in hits:
            all_rows.append(dict(
                protein=pid, isoform=h["iso"], dual=h["dual"],
                win_start=h["win_start"], win_end=h["win_end"],
                core=h["core"], window=h["window"], pssm_score=h["pssm_score"],
                canonical=h["canonical"], amphipathic=h["amphipathic"],
                pI=h["pI"], helix_approx=h["helix_approx"]))

    all_rows.sort(key=lambda r: -r["pssm_score"])

    # --- Write output ---
    fields = ["protein", "isoform", "dual", "win_start", "win_end", "core", "window",
              "pssm_score", "canonical", "amphipathic", "pI", "helix_approx"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    hit_proteins = sorted(set(r["protein"] for r in all_rows))
    no_hit = sorted(set(pid for pid, _ in entries) - set(hit_proteins))

    # --- Print summary ---
    print("=" * 60)
    print("  RESULTS")
    print("=" * 60)
    print(f"\n  Proteins screened:  {len(entries)}")
    print(f"  Proteins with hit:  {len(hit_proteins)}")
    print(f"  Proteins no hit:    {len(no_hit)}")
    print(f"  Total motif hits:   {len(all_rows)}")

    if hit_proteins:
        print(f"\n  Proteins WITH AKAP motif(s):")
        for pid in hit_proteins:
            hits_for = [r for r in all_rows if r["protein"] == pid]
            isos = "+".join(sorted(set(r["isoform"] for r in hits_for)))
            dual = any(r["dual"] for r in hits_for)
            best = max(r["pssm_score"] for r in hits_for)
            tag = " (dual-specific)" if dual else ""
            print(f"    {pid:30s}  {isos:6s}  best_score={best:6.2f}{tag}")

    if no_hit:
        print(f"\n  Proteins WITHOUT AKAP motif:")
        for pid in no_hit:
            print(f"    {pid}")

    print(f"\n  Output saved to: {args.out}")
    if args.keep_fasta or args.fasta_out:
        print(f"  FASTA saved to:  {fasta_path}")
    else:
        os.unlink(fasta_path)

    print()


if __name__ == "__main__":
    main()
