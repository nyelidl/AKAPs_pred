#!/usr/bin/env python3
"""
build_candidate_list.py
Defines the biological hard-negative candidate set (accessions + context only,
NO sequences) and writes:
  - hard_negative_candidate_list_biological_context.csv
  - hard_negative_candidate_list_biological_context.md
Every row carries full biological-context annotation. `verification_status`
flags accessions/regions that MUST be confirmed before fetching sequence.
Ambiguous-AKAP-status candidates are marked priority=flagged_ambiguous and are
NOT to be used silently as negatives.
"""
import csv, os

FIELDS = ["id","class","protein_or_peptide_name","organism","accession","database",
          "pdb_id","sequence_type","known_region","cellular_location",
          "biological_function","known_binding_partner","known_binding_surface",
          "reason_akap_like","reason_not_akap","biological_context_evidence",
          "structural_context_evidence","expected_failure_mode","intended_use",
          "priority","verification_status","notes"]

def C(**kw):
    return {f: kw.get(f, "") for f in FIELDS}

candidates = [
# ───────────────── Class 1 — Coiled-coil amphipathic helices ─────────────────
C(id="CC01", **{"class":"coiled_coil"}, protein_or_peptide_name="GCN4 leucine zipper",
  organism="Saccharomyces cerevisiae", accession="P03069", database="UniProt", pdb_id="2ZTA",
  sequence_type="full_protein", known_region="249-281", cellular_location="nucleus",
  biological_function="bZIP transcription factor homodimerization",
  known_binding_partner="GCN4 leucine zipper (homodimer)",
  known_binding_surface="helix-helix coiled-coil interface (heptad a/d)",
  reason_akap_like="continuous hydrophobic seam, amphipathic, can score on hydrophobic stripe",
  reason_not_akap="hydrophobic face packs against a partner helix (knobs-into-holes), not the PKA D/D groove; 3.5 res/turn heptad",
  biological_context_evidence="canonical bZIP dimerization domain (O'Shea 1991)",
  structural_context_evidence="PDB 2ZTA parallel dimeric coiled coil",
  expected_failure_mode="hydrophobic stripe may give PSSM>=12; ML may follow PSSM",
  intended_use="training_candidate+stress_test", priority="high",
  verification_status="verify_region", notes="prototypical coiled coil; clean label"),
C(id="CC02", **{"class":"coiled_coil"}, protein_or_peptide_name="Tropomyosin alpha-1 chain",
  organism="Homo sapiens", accession="P09493", database="UniProt", pdb_id="1C1G",
  sequence_type="full_protein", known_region="full coiled-coil (use sliding window)",
  cellular_location="cytoskeleton",
  biological_function="actin filament stabilization; coiled-coil dimer",
  known_binding_partner="tropomyosin partner chain / actin",
  known_binding_surface="helix-helix coiled-coil interface",
  reason_akap_like="long amphipathic coiled coil with strong hydrophobic seam",
  reason_not_akap="obligate dimerization seam; not a PKA anchor",
  biological_context_evidence="muscle thin-filament regulator",
  structural_context_evidence="PDB 1C1G",
  expected_failure_mode="multiple windows may breach danger zone along the rod",
  intended_use="training_candidate+stress_test", priority="high",
  verification_status="verify_accession_isoform", notes="very long; many windows"),
C(id="CC03", **{"class":"coiled_coil"}, protein_or_peptide_name="Vimentin (coil 2B rod)",
  organism="Homo sapiens", accession="P08670", database="UniProt", pdb_id="3UF1",
  sequence_type="full_protein", known_region="coil 2B (~262-334)", cellular_location="cytoskeleton",
  biological_function="intermediate filament assembly",
  known_binding_partner="vimentin (dimer/tetramer)",
  known_binding_surface="helix-helix coiled-coil interface",
  reason_akap_like="amphipathic coiled-coil rod",
  reason_not_akap="IF assembly interface, not PKA D/D",
  biological_context_evidence="canonical type-III IF",
  structural_context_evidence="PDB 3UF1 coil2 fragment",
  expected_failure_mode="danger-zone entry along rod", intended_use="stress_test",
  priority="medium", verification_status="verify_region", notes=""),
C(id="CC04", **{"class":"coiled_coil"}, protein_or_peptide_name="SNAP-25 SNARE motif",
  organism="Homo sapiens", accession="P60880", database="UniProt", pdb_id="1SFC",
  sequence_type="full_protein", known_region="SN1/SN2 SNARE helices", cellular_location="membrane/vesicle",
  biological_function="SNARE complex membrane fusion",
  known_binding_partner="syntaxin-1 / VAMP2 (SNARE 4-helix bundle)",
  known_binding_surface="four-helix bundle interface",
  reason_akap_like="amphipathic helices with hydrophobic interface",
  reason_not_akap="assembles into SNARE bundle; not PKA anchor",
  biological_context_evidence="neuronal exocytosis SNARE",
  structural_context_evidence="PDB 1SFC", expected_failure_mode="bundle-facing hydrophobic face may score",
  intended_use="stress_test", priority="medium", verification_status="verify_region", notes=""),

# ─────────────── Class 2 — Membrane-binding amphipathic helices ───────────────
C(id="MB01", **{"class":"membrane_amphipathic"}, protein_or_peptide_name="ArfGAP1 ALPS1 motif",
  organism="Homo sapiens", accession="Q8N6T3", database="UniProt", pdb_id="",
  sequence_type="full_protein", known_region="ALPS1 ~192-257 (lipid packing sensor)",
  cellular_location="Golgi membrane / cytosol",
  biological_function="curvature-sensing membrane binding (ArfGAP activity)",
  known_binding_partner="curved lipid membrane",
  known_binding_surface="membrane-facing hydrophobic face",
  reason_akap_like="strong amphipathic helix, high hydrophobic moment",
  reason_not_akap="inserts bulky face into lipid packing defects; charge-poor polar face; not PKA groove",
  biological_context_evidence="ALPS curvature sensor (Bigay/Antonny)",
  structural_context_evidence="folds on membrane; disordered in solution",
  expected_failure_mode="HIGH likelihood of danger-zone entry; key test class",
  intended_use="training_candidate+stress_test", priority="high",
  verification_status="verify_accession_region", notes="charge-poor polar face is the discriminator"),
C(id="MB02", **{"class":"membrane_amphipathic"}, protein_or_peptide_name="Alpha-synuclein N-terminal repeats",
  organism="Homo sapiens", accession="P37840", database="UniProt", pdb_id="1XQ8",
  sequence_type="full_protein", known_region="1-95 (KTKEGV imperfect 11-mers)",
  cellular_location="cytosol / presynaptic membrane",
  biological_function="synaptic vesicle membrane binding",
  known_binding_partner="acidic phospholipid membrane",
  known_binding_surface="membrane-facing amphipathic face",
  reason_akap_like="long amphipathic helix on membrane",
  reason_not_akap="lipid-binding; 11-mer periodicity; lysine-rich interface; not PKA D/D",
  biological_context_evidence="micelle-bound helix (Ulmer 2005)",
  structural_context_evidence="PDB 1XQ8 (micelle-bound)",
  expected_failure_mode="danger-zone entry likely across repeats",
  intended_use="training_candidate+stress_test", priority="high",
  verification_status="verify_region", notes=""),
C(id="MB03", **{"class":"membrane_amphipathic"}, protein_or_peptide_name="Epsin-1 ENTH N-terminal helix-0",
  organism="Homo sapiens", accession="Q9Y6I3", database="UniProt", pdb_id="1H0A",
  sequence_type="full_protein", known_region="helix-0 (~1-15)", cellular_location="plasma membrane / cytosol",
  biological_function="PtdIns(4,5)P2-induced membrane curvature",
  known_binding_partner="PIP2 membrane",
  known_binding_surface="membrane-inserted hydrophobic face",
  reason_akap_like="amphipathic helix that folds on membrane",
  reason_not_akap="short curvature-driving helix; lipid target; not PKA",
  biological_context_evidence="ENTH H0 insertion (Ford 2002)",
  structural_context_evidence="PDB 1H0A", expected_failure_mode="short; may not reach 5-turn support",
  intended_use="stress_test", priority="medium", verification_status="verify_accession_region",
  notes="tests 5-turn/length discriminator"),
C(id="MB04", **{"class":"membrane_amphipathic"}, protein_or_peptide_name="Cathelicidin LL-37",
  organism="Homo sapiens", accession="P49913", database="UniProt", pdb_id="2K6O",
  sequence_type="peptide_region", known_region="mature 134-170", cellular_location="secreted",
  biological_function="antimicrobial membrane disruption",
  known_binding_partner="microbial membrane",
  known_binding_surface="membrane-facing hydrophobic face",
  reason_akap_like="cationic amphipathic helix, high moment",
  reason_not_akap="membrane-lytic; strongly cationic polar face; not PKA groove",
  biological_context_evidence="human cathelicidin AMP",
  structural_context_evidence="PDB 2K6O", expected_failure_mode="cationic face; PSSM may reject",
  intended_use="stress_test", priority="medium", verification_status="verify_region",
  notes="peptide-derived region; cationic polar face discriminator"),
C(id="MB05", **{"class":"membrane_amphipathic"}, protein_or_peptide_name="Melittin",
  organism="Apis mellifera", accession="P01501", database="UniProt", pdb_id="2MLT",
  sequence_type="peptide_only", known_region="1-26 (mature)", cellular_location="secreted (venom)",
  biological_function="membrane-lytic peptide",
  known_binding_partner="lipid membrane", known_binding_surface="membrane-facing hydrophobic face",
  reason_akap_like="classic amphipathic helix, very high moment",
  reason_not_akap="pore-forming lytic peptide; not a protein-groove anchor",
  biological_context_evidence="bee venom peptide",
  structural_context_evidence="PDB 2MLT", expected_failure_mode="peptide-only; lacks protein context",
  intended_use="stress_test_only", priority="low", verification_status="verify_sequence",
  notes="PEPTIDE-ONLY: no full-protein background; analyze separately"),

# ────────────── Class 3 — Mitochondrial targeting presequences ──────────────
C(id="MT01", **{"class":"mito_targeting"}, protein_or_peptide_name="ALDH2 presequence",
  organism="Homo sapiens", accession="P05091", database="UniProt", pdb_id="",
  sequence_type="full_protein", known_region="N-terminal presequence ~1-17",
  cellular_location="mitochondria (matrix)",
  biological_function="mitochondrial import targeting (cleaved)",
  known_binding_partner="TOM20 import receptor",
  known_binding_surface="TOM20 hydrophobic groove (phi-chi-chi-phi-phi)",
  reason_akap_like="amphipathic with hydrophobic + charged face",
  reason_not_akap="Arg-rich, acidic-depleted targeting peptide; cleaved; recognized by TOM/TIM not PKA",
  biological_context_evidence="classic cleavable MTS",
  structural_context_evidence="TOM20-presequence complexes (Saitoh 2007)",
  expected_failure_mode="likely PSSM<12 (different anchors) -> PSSM filters",
  intended_use="stress_test", priority="medium", verification_status="verify_region",
  notes="if PSSM<12, report PSSM already filters this class"),
C(id="MT02", **{"class":"mito_targeting"}, protein_or_peptide_name="Ornithine carbamoyltransferase presequence",
  organism="Homo sapiens", accession="P00480", database="UniProt", pdb_id="",
  sequence_type="full_protein", known_region="N-terminal presequence ~1-32",
  cellular_location="mitochondria (matrix)", biological_function="import targeting (cleaved)",
  known_binding_partner="TOM20", known_binding_surface="TOM20 groove",
  reason_akap_like="amphipathic targeting helix",
  reason_not_akap="positive-inside, acidic-poor; TOM/TIM substrate",
  biological_context_evidence="canonical OTC presequence",
  structural_context_evidence="", expected_failure_mode="likely PSSM<12",
  intended_use="stress_test", priority="low", verification_status="verify_region", notes=""),
C(id="MT03", **{"class":"mito_targeting"}, protein_or_peptide_name="Cytochrome c oxidase subunit 4 isoform 1 presequence",
  organism="Homo sapiens", accession="P13073", database="UniProt", pdb_id="",
  sequence_type="full_protein", known_region="N-terminal presequence ~1-22",
  cellular_location="mitochondria", biological_function="import targeting (cleaved)",
  known_binding_partner="TOM20", known_binding_surface="TOM20 groove",
  reason_akap_like="amphipathic presequence helix",
  reason_not_akap="Arg-rich acidic-poor targeting peptide",
  biological_context_evidence="model MTS (COX4)", structural_context_evidence="",
  expected_failure_mode="likely PSSM<12", intended_use="stress_test", priority="low",
  verification_status="verify_region", notes=""),

# ──────────── Class 4 — Lipid-binding / apolipoprotein helices ────────────
C(id="LP01", **{"class":"apolipoprotein"}, protein_or_peptide_name="Apolipoprotein A-I",
  organism="Homo sapiens", accession="P02647", database="UniProt", pdb_id="1AV1",
  sequence_type="full_protein", known_region="class-A 11/22-mer repeats (mature)",
  cellular_location="secreted (HDL)", biological_function="lipid binding / HDL scaffold",
  known_binding_partner="phospholipid / lipoprotein surface",
  known_binding_surface="lipid-facing hydrophobic face",
  reason_akap_like="textbook class-A amphipathic helix, sharp face separation, high moment",
  reason_not_akap="lipid binder; interfacial-basic + center-acidic snorkel; tandem repeats; not PKA",
  biological_context_evidence="major HDL apolipoprotein",
  structural_context_evidence="PDB 1AV1 (lipid-free belt)",
  expected_failure_mode="HIGH danger-zone entry; central acidic cluster may trigger negdet->unlikely",
  intended_use="training_candidate+stress_test", priority="high",
  verification_status="verify_region", notes="also tests whether negdet rule earns its keep"),
C(id="LP02", **{"class":"apolipoprotein"}, protein_or_peptide_name="Apolipoprotein E (C-terminal lipid-binding)",
  organism="Homo sapiens", accession="P02649", database="UniProt", pdb_id="1B68",
  sequence_type="full_protein", known_region="C-terminal lipid-binding helices ~244-272",
  cellular_location="secreted", biological_function="lipoprotein lipid binding",
  known_binding_partner="lipoprotein lipid surface", known_binding_surface="lipid-facing face",
  reason_akap_like="amphipathic lipid-binding helices",
  reason_not_akap="lipid target; not PKA D/D", biological_context_evidence="ApoE lipid binding",
  structural_context_evidence="PDB 1B68 region", expected_failure_mode="danger-zone entry plausible",
  intended_use="stress_test", priority="medium", verification_status="verify_region", notes=""),
C(id="LP03", **{"class":"apolipoprotein"}, protein_or_peptide_name="Perilipin-1 (class-A repeats)",
  organism="Homo sapiens", accession="O60240", database="UniProt", pdb_id="",
  sequence_type="full_protein", known_region="11-mer amphipathic repeat region",
  cellular_location="lipid droplet", biological_function="lipid droplet coat",
  known_binding_partner="lipid droplet monolayer", known_binding_surface="lipid-facing face",
  reason_akap_like="class-A amphipathic 11-mer repeats",
  reason_not_akap="lipid-droplet binder; not PKA (note: PLIN1 is regulated by PKA but anchoring helix is lipid-facing)",
  biological_context_evidence="PAT-family LD protein",
  structural_context_evidence="", expected_failure_mode="danger-zone entry plausible",
  intended_use="stress_test", priority="medium", verification_status="verify_region",
  notes="PLIN1 is a PKA SUBSTRATE - ensure region chosen is the lipid-binding helix, not a PKA site"),

# ───────────── Class 5 — Generic helix-in-groove PPI motifs ─────────────
C(id="PG01", **{"class":"ppi_helix_in_groove"}, protein_or_peptide_name="BAD BH3 helix",
  organism="Homo sapiens", accession="Q92934", database="UniProt", pdb_id="1G5J",
  sequence_type="full_protein", known_region="BH3 ~103-127", cellular_location="cytosol/mitochondria",
  biological_function="apoptosis regulation (pro-apoptotic BH3-only)",
  known_binding_partner="Bcl-XL / Bcl-2 hydrophobic groove",
  known_binding_surface="protein-groove interface",
  reason_akap_like="short amphipathic helix docking into a hydrophobic protein groove (mechanistically closest)",
  reason_not_akap="BH3 anchor spacing (Lxxxx(D/phi)) and Bcl-2 partner differ from AKAP D/D",
  biological_context_evidence="BH3-in-groove binding",
  structural_context_evidence="PDB 1G5J", expected_failure_mode="anchor spacing differs; possible PSSM<12",
  intended_use="stress_test", priority="high", verification_status="verify_region",
  notes="most diagnostic if it reaches danger zone"),
C(id="PG02", **{"class":"ppi_helix_in_groove"}, protein_or_peptide_name="BIM BH3 helix",
  organism="Homo sapiens", accession="O43521", database="UniProt", pdb_id="1PQ1",
  sequence_type="full_protein", known_region="BH3 ~141-166", cellular_location="cytosol/mitochondria",
  biological_function="apoptosis (BH3-only)", known_binding_partner="Bcl-2 family groove",
  known_binding_surface="protein-groove interface", reason_akap_like="helix-in-groove amphipathic",
  reason_not_akap="Bcl-2 partner; BH3 motif spacing", biological_context_evidence="BH3 binding",
  structural_context_evidence="PDB 1PQ1", expected_failure_mode="anchor spacing differs",
  intended_use="stress_test", priority="medium", verification_status="verify_region", notes=""),
C(id="PG03", **{"class":"ppi_helix_in_groove"}, protein_or_peptide_name="p53 transactivation helix",
  organism="Homo sapiens", accession="P04637", database="UniProt", pdb_id="1YCR",
  sequence_type="full_protein", known_region="TAD helix ~17-29", cellular_location="nucleus",
  biological_function="transcription; MDM2 regulation",
  known_binding_partner="MDM2 hydrophobic cleft", known_binding_surface="protein-groove interface",
  reason_akap_like="short helix in hydrophobic groove (F19/W23/L26 anchors)",
  reason_not_akap="MDM2 partner; 3-residue anchor spacing; very short (~3 turns)",
  biological_context_evidence="p53-MDM2 (Kussie 1996)", structural_context_evidence="PDB 1YCR",
  expected_failure_mode="too short for 5-turn AKAP support; likely filtered",
  intended_use="stress_test", priority="high", verification_status="verify_region",
  notes="tests 5-turn support discriminator"),
C(id="PG04", **{"class":"ppi_helix_in_groove"}, protein_or_peptide_name="NCOA1/SRC-1 LxxLL NR box",
  organism="Homo sapiens", accession="Q15788", database="UniProt", pdb_id="1GWQ",
  sequence_type="full_protein", known_region="NR box LxxLL motif", cellular_location="nucleus",
  biological_function="nuclear receptor coactivation",
  known_binding_partner="nuclear receptor LBD AF-2 groove",
  known_binding_surface="protein-groove interface",
  reason_akap_like="short amphipathic helix in a hydrophobic groove",
  reason_not_akap="LxxLL motif; NR-LBD partner; short", biological_context_evidence="NR box binding",
  structural_context_evidence="PDB 1GWQ", expected_failure_mode="short motif; likely filtered",
  intended_use="stress_test", priority="medium", verification_status="verify_region", notes=""),

# ───────────── Class 6 — DDIP-like / DPY30-binding helices ─────────────
C(id="DD01", **{"class":"ddip_like"}, protein_or_peptide_name="DPY30",
  organism="Homo sapiens", accession="Q9C005", database="UniProt", pdb_id="3G36",
  sequence_type="full_protein", known_region="dimerization/binding region", cellular_location="nucleus",
  biological_function="MLL/COMPASS assembly; binds D/D-like fold",
  known_binding_partner="ASH2L (Sdc-binding) / D/D-like domain",
  known_binding_surface="D/D-like domain interface (NOT PKA R-subunit)",
  reason_akap_like="engages a dimerization/docking (D/D-like) fold via helix",
  reason_not_akap="partner is DPY30/ASH2L system, not PKA RI/RII; architecture differs",
  biological_context_evidence="COMPASS subunit",
  structural_context_evidence="PDB 3G36 region",
  expected_failure_mode="D/D-like but non-PKA; tests DDIP discrimination",
  intended_use="stress_test", priority="high",
  verification_status="verify_region_and_AKAP_status",
  notes="FLAGGED: confirm it is experimentally NOT a PKA anchor before training use"),
C(id="DD02", **{"class":"ddip_like"}, protein_or_peptide_name="ASH2L (DPY30-binding helix)",
  organism="Homo sapiens", accession="Q9UBL3", database="UniProt", pdb_id="",
  sequence_type="full_protein", known_region="DPY30-binding region", cellular_location="nucleus",
  biological_function="COMPASS H3K4 methyltransferase complex",
  known_binding_partner="DPY30 D/D-like domain", known_binding_surface="D/D-like interface",
  reason_akap_like="helix binding a D/D-like fold", reason_not_akap="non-PKA D/D-like partner",
  biological_context_evidence="ASH2L-DPY30 interaction",
  structural_context_evidence="", expected_failure_mode="DDIP-like discrimination test",
  intended_use="stress_test", priority="medium",
  verification_status="verify_region_and_AKAP_status",
  notes="FLAGGED ambiguous: confirm non-PKA-anchoring before any training use"),
]

OUT = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(OUT, "hard_negative_candidate_list_biological_context.csv")
with open(csv_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader()
    for c in candidates:
        w.writerow(c)

# human-readable markdown grouped by class
md_path = os.path.join(OUT, "hard_negative_candidate_list_biological_context.md")
by_class = {}
for c in candidates:
    by_class.setdefault(c["class"], []).append(c)
CLASS_TITLE = {
    "coiled_coil":"Class 1 — Coiled-coil amphipathic helices",
    "membrane_amphipathic":"Class 2 — Membrane-binding amphipathic helices",
    "mito_targeting":"Class 3 — Mitochondrial targeting presequences",
    "apolipoprotein":"Class 4 — Lipid-binding / apolipoprotein helices",
    "ppi_helix_in_groove":"Class 5 — Generic helix-in-groove PPI motifs",
    "ddip_like":"Class 6 — DDIP-like / DPY30-binding helices",
}
with open(md_path, "w") as fh:
    fh.write("# Hard-negative candidate list (biological context)\n\n")
    fh.write("Accessions and regions only — **no sequences**. All entries are "
             "`label=non_AKAP`. `verification_status` MUST be cleared before fetching "
             "sequence. Entries flagged `verify_..._AKAP_status` are ambiguous and must "
             "NOT be used as training negatives until their non-PKA-anchoring status is "
             "confirmed.\n\n")
    fh.write(f"Total candidates: {len(candidates)} across {len(by_class)} classes.\n\n")
    for cls in CLASS_TITLE:
        rows = by_class.get(cls, [])
        if not rows: continue
        fh.write(f"## {CLASS_TITLE[cls]}\n\n")
        for c in rows:
            fh.write(f"### {c['id']} — {c['protein_or_peptide_name']} ({c['organism']})\n")
            fh.write(f"- accession: **{c['accession']}** ({c['database']}); PDB: {c['pdb_id'] or '—'}; "
                     f"region: {c['known_region']}; type: {c['sequence_type']}\n")
            fh.write(f"- location: {c['cellular_location']} | function: {c['biological_function']}\n")
            fh.write(f"- binding partner: {c['known_binding_partner']} | surface: {c['known_binding_surface']}\n")
            fh.write(f"- AKAP-like because: {c['reason_akap_like']}\n")
            fh.write(f"- NOT AKAP because: {c['reason_not_akap']}\n")
            fh.write(f"- context evidence: {c['biological_context_evidence']}; "
                     f"structural: {c['structural_context_evidence'] or '—'}\n")
            fh.write(f"- expected failure mode: {c['expected_failure_mode']}\n")
            fh.write(f"- intended use: **{c['intended_use']}** | priority: **{c['priority']}** | "
                     f"verification: {c['verification_status']}\n")
            if c['notes']:
                fh.write(f"- notes: {c['notes']}\n")
            fh.write("\n")

print(f"wrote {csv_path}")
print(f"wrote {md_path}")
print(f"candidates: {len(candidates)} | classes: {sorted(by_class)}")
for cls in by_class:
    print(f"  {cls}: {len(by_class[cls])}")
