#!/usr/bin/env python3
"""
hard_negative_qc.py
QC for the biological hard-negative benchmark. Runs on the metadata/candidate CSV
and (optionally) the verified FASTA. Enforces the curation rules BEFORE screening:

  1. FASTA IDs match metadata IDs (both directions), if a FASTA is given.
  2. Metadata has the required biological-context fields, non-empty.
  3. No sequence accepted without accession OR source evidence.
  4. Peptide-only examples are labelled separately (sequence_type).
  5. Full-protein examples flagged as preferred for deployment realism.
  6. set_type groups (positive_control / easy_negative_control /
     synthetic_mechanistic_decoy / biological_hard_negative) are not mixed.
  7. Each biological hard negative has BOTH reason_akap_like and reason_not_akap.
  8. Ambiguous-AKAP-status candidates are flagged, not silently used as
     training negatives (verification_status mentioning AKAP_status must not
     have intended_use containing 'training').

Exit code 0 = pass (warnings allowed), 1 = hard failures present.

Usage:
  python3 hard_negative_qc.py --meta hard_negative_candidate_list_biological_context.csv
  python3 hard_negative_qc.py --meta metadata.csv --fasta hard_negative_amphipathic_set.fasta
"""
import argparse, csv, os, sys

REQUIRED_CONTEXT = ["class","protein_or_peptide_name","organism","accession","database",
                    "sequence_type","known_region","cellular_location","biological_function",
                    "known_binding_partner","known_binding_surface","reason_akap_like",
                    "reason_not_akap","biological_context_evidence","intended_use","priority"]
VALID_SET_TYPES = {"biological_hard_negative","positive_control",
                   "easy_negative_control","synthetic_mechanistic_decoy"}

def read_meta(path):
    with open(path) as fh:
        return [r for r in csv.DictReader(fh) if any(v.strip() for v in r.values())]

def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0].split("|")[0])
    return ids

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--fasta", default=None)
    a = ap.parse_args()

    errors, warnings = [], []
    rows = read_meta(a.meta)
    if not rows:
        print("FAIL: metadata empty"); sys.exit(1)

    rows = [r for r in rows if r.get("id","").strip() != "EXAMPLE_DELETE_ME"]
    ids = [r["id"].strip() for r in rows]

    # 0. unique IDs
    dups = {i for i in ids if ids.count(i) > 1}
    if dups: errors.append(f"[unique-id] duplicate IDs: {sorted(dups)}")

    has_set_type = "set_type" in rows[0]
    for r in rows:
        rid = r.get("id","?")
        st = (r.get("set_type","") or "biological_hard_negative").strip()
        is_bio_hn = (st == "biological_hard_negative")

        # 6. set_type validity
        if has_set_type and st and st not in VALID_SET_TYPES:
            errors.append(f"[{rid}] invalid set_type '{st}'")

        # 2 + 7. required context fields (biological hard negatives)
        if is_bio_hn:
            for f in REQUIRED_CONTEXT:
                if not (r.get(f,"") or "").strip():
                    errors.append(f"[{rid}] missing required context field '{f}'")

        # 3. accession or source evidence
        if not (r.get("accession","").strip() or r.get("biological_context_evidence","").strip()
                or r.get("structural_context_evidence","").strip()):
            errors.append(f"[{rid}] no accession AND no source evidence")

        # 1b. label must be non_AKAP for biological hard negatives
        if is_bio_hn and "label" in r and r.get("label","").strip() not in ("", "non_AKAP"):
            errors.append(f"[{rid}] biological hard negative must be label=non_AKAP, got '{r['label']}'")

        # 4 + 5. peptide-only vs full-protein
        stype = r.get("sequence_type","").strip()
        if stype == "peptide_only":
            warnings.append(f"[{rid}] peptide-only: analyse separately (no full-protein background)")
        elif stype == "full_protein":
            pass  # preferred
        elif is_bio_hn and stype == "":
            errors.append(f"[{rid}] sequence_type missing (need full_protein / peptide_region / peptide_only)")

        # 8. ambiguous AKAP status must not be a training negative
        vstat = r.get("verification_status","").lower()
        iuse  = r.get("intended_use","").lower()
        if "akap_status" in vstat and "training" in iuse:
            errors.append(f"[{rid}] ambiguous AKAP status but intended_use includes training "
                          f"-> must be stress_test/holdout only until status confirmed")
        if "akap_status" in vstat:
            warnings.append(f"[{rid}] AKAP status unconfirmed -> flagged, not silently used as negative")

        # verification gate before fetch/use
        if r.get("verification_status","").strip() not in ("", "verified") and a.fasta:
            warnings.append(f"[{rid}] verification_status='{r['verification_status']}' but a FASTA was supplied "
                            f"-> confirm before trusting this row")

    # 1. FASTA <-> metadata ID parity
    if a.fasta:
        if not os.path.exists(a.fasta):
            errors.append(f"[fasta] not found: {a.fasta}")
        else:
            fids = read_fasta_ids(a.fasta)
            fset, mset = set(fids), set(ids)
            for x in fset - mset: errors.append(f"[fasta] id '{x}' has no metadata row")
            for x in mset - fset: warnings.append(f"[meta] id '{x}' has no FASTA entry (not yet fetched?)")
            d = {i for i in fids if fids.count(i) > 1}
            if d: errors.append(f"[fasta] duplicate FASTA ids: {sorted(d)}")

    # set_type mixing note
    if has_set_type:
        from collections import Counter
        c = Counter((r.get("set_type","") or "biological_hard_negative").strip() for r in rows)
        print("set_type composition:", dict(c))
        if not has_set_type:
            warnings.append("[set_type] column absent -> cannot guarantee controls are separated")

    print(f"\n{len(rows)} rows checked.")
    if warnings:
        print(f"\n{len(warnings)} WARNING(S):")
        for w in warnings: print("  ! "+w)
    if errors:
        print(f"\n{len(errors)} ERROR(S):")
        for e in errors: print("  x "+e)
        print("\nQC RESULT: FAIL"); sys.exit(1)
    print("\nQC RESULT: PASS"); sys.exit(0)

if __name__ == "__main__":
    main()
