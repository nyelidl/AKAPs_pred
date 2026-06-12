# Hard-negative candidate list (biological context)

Accessions and regions only — **no sequences**. All entries are `label=non_AKAP`. `verification_status` MUST be cleared before fetching sequence. Entries flagged `verify_..._AKAP_status` are ambiguous and must NOT be used as training negatives until their non-PKA-anchoring status is confirmed.

Total candidates: 25 across 7 classes.

## Class 1 — Coiled-coil amphipathic helices

### CC01 — GCN4 leucine zipper (Saccharomyces cerevisiae)
- accession: **P03069** (UniProt); PDB: 2ZTA; region: 249-281; type: full_protein
- location: nucleus | function: bZIP transcription factor homodimerization
- binding partner: GCN4 leucine zipper (homodimer) | surface: helix-helix coiled-coil interface (heptad a/d)
- AKAP-like because: continuous hydrophobic seam, amphipathic, can score on hydrophobic stripe
- NOT AKAP because: hydrophobic face packs against a partner helix (knobs-into-holes), not the PKA D/D groove; 3.5 res/turn heptad
- context evidence: canonical bZIP dimerization domain (O'Shea 1991); structural: PDB 2ZTA parallel dimeric coiled coil
- expected failure mode: hydrophobic stripe may give PSSM>=12; ML may follow PSSM
- intended use: **training_candidate+stress_test** | priority: **high** | verification: verify_region
- notes: prototypical coiled coil; clean label

### CC02 — Tropomyosin alpha-1 chain (Homo sapiens)
- accession: **P09493** (UniProt); PDB: 1C1G; region: full coiled-coil (use sliding window); type: full_protein
- location: cytoskeleton | function: actin filament stabilization; coiled-coil dimer
- binding partner: tropomyosin partner chain / actin | surface: helix-helix coiled-coil interface
- AKAP-like because: long amphipathic coiled coil with strong hydrophobic seam
- NOT AKAP because: obligate dimerization seam; not a PKA anchor
- context evidence: muscle thin-filament regulator; structural: PDB 1C1G
- expected failure mode: multiple windows may breach danger zone along the rod
- intended use: **training_candidate+stress_test** | priority: **high** | verification: verify_accession_isoform
- notes: very long; many windows

### CC03 — Vimentin (coil 2B rod) (Homo sapiens)
- accession: **P08670** (UniProt); PDB: 3UF1; region: coil 2B (~262-334); type: full_protein
- location: cytoskeleton | function: intermediate filament assembly
- binding partner: vimentin (dimer/tetramer) | surface: helix-helix coiled-coil interface
- AKAP-like because: amphipathic coiled-coil rod
- NOT AKAP because: IF assembly interface, not PKA D/D
- context evidence: canonical type-III IF; structural: PDB 3UF1 coil2 fragment
- expected failure mode: danger-zone entry along rod
- intended use: **stress_test** | priority: **medium** | verification: verify_region

### CC04 — SNAP-25 SNARE motif (Homo sapiens)
- accession: **P60880** (UniProt); PDB: 1SFC; region: SN1/SN2 SNARE helices; type: full_protein
- location: membrane/vesicle | function: SNARE complex membrane fusion
- binding partner: syntaxin-1 / VAMP2 (SNARE 4-helix bundle) | surface: four-helix bundle interface
- AKAP-like because: amphipathic helices with hydrophobic interface
- NOT AKAP because: assembles into SNARE bundle; not PKA anchor
- context evidence: neuronal exocytosis SNARE; structural: PDB 1SFC
- expected failure mode: bundle-facing hydrophobic face may score
- intended use: **stress_test** | priority: **medium** | verification: verify_region

## Class 2 — Membrane-binding amphipathic helices

### MB01 — ArfGAP1 ALPS1 motif (Homo sapiens)
- accession: **Q8N6T3** (UniProt); PDB: —; region: ALPS1 ~192-257 (lipid packing sensor); type: full_protein
- location: Golgi membrane / cytosol | function: curvature-sensing membrane binding (ArfGAP activity)
- binding partner: curved lipid membrane | surface: membrane-facing hydrophobic face
- AKAP-like because: strong amphipathic helix, high hydrophobic moment
- NOT AKAP because: inserts bulky face into lipid packing defects; charge-poor polar face; not PKA groove
- context evidence: ALPS curvature sensor (Bigay/Antonny); structural: folds on membrane; disordered in solution
- expected failure mode: HIGH likelihood of danger-zone entry; key test class
- intended use: **training_candidate+stress_test** | priority: **high** | verification: verify_accession_region
- notes: charge-poor polar face is the discriminator

### MB02 — Alpha-synuclein N-terminal repeats (Homo sapiens)
- accession: **P37840** (UniProt); PDB: 1XQ8; region: 1-95 (KTKEGV imperfect 11-mers); type: full_protein
- location: cytosol / presynaptic membrane | function: synaptic vesicle membrane binding
- binding partner: acidic phospholipid membrane | surface: membrane-facing amphipathic face
- AKAP-like because: long amphipathic helix on membrane
- NOT AKAP because: lipid-binding; 11-mer periodicity; lysine-rich interface; not PKA D/D
- context evidence: micelle-bound helix (Ulmer 2005); structural: PDB 1XQ8 (micelle-bound)
- expected failure mode: danger-zone entry likely across repeats
- intended use: **training_candidate+stress_test** | priority: **high** | verification: verify_region

### MB03 — Epsin-1 ENTH N-terminal helix-0 (Homo sapiens)
- accession: **Q9Y6I3** (UniProt); PDB: 1H0A; region: helix-0 (~1-15); type: full_protein
- location: plasma membrane / cytosol | function: PtdIns(4,5)P2-induced membrane curvature
- binding partner: PIP2 membrane | surface: membrane-inserted hydrophobic face
- AKAP-like because: amphipathic helix that folds on membrane
- NOT AKAP because: short curvature-driving helix; lipid target; not PKA
- context evidence: ENTH H0 insertion (Ford 2002); structural: PDB 1H0A
- expected failure mode: short; may not reach 5-turn support
- intended use: **stress_test** | priority: **medium** | verification: verify_accession_region
- notes: tests 5-turn/length discriminator

### MB04 — Cathelicidin LL-37 (Homo sapiens)
- accession: **P49913** (UniProt); PDB: 2K6O; region: mature 134-170; type: peptide_region
- location: secreted | function: antimicrobial membrane disruption
- binding partner: microbial membrane | surface: membrane-facing hydrophobic face
- AKAP-like because: cationic amphipathic helix, high moment
- NOT AKAP because: membrane-lytic; strongly cationic polar face; not PKA groove
- context evidence: human cathelicidin AMP; structural: PDB 2K6O
- expected failure mode: cationic face; PSSM may reject
- intended use: **stress_test** | priority: **medium** | verification: verify_region
- notes: peptide-derived region; cationic polar face discriminator

### MB05 — Melittin (Apis mellifera)
- accession: **P01501** (UniProt); PDB: 2MLT; region: 1-26 (mature); type: peptide_only
- location: secreted (venom) | function: membrane-lytic peptide
- binding partner: lipid membrane | surface: membrane-facing hydrophobic face
- AKAP-like because: classic amphipathic helix, very high moment
- NOT AKAP because: pore-forming lytic peptide; not a protein-groove anchor
- context evidence: bee venom peptide; structural: PDB 2MLT
- expected failure mode: peptide-only; lacks protein context
- intended use: **stress_test_only** | priority: **low** | verification: verify_sequence
- notes: PEPTIDE-ONLY: no full-protein background; analyze separately

## Class 3 — Mitochondrial targeting presequences

### MT01 — ALDH2 presequence (Homo sapiens)
- accession: **P05091** (UniProt); PDB: —; region: N-terminal presequence ~1-17; type: full_protein
- location: mitochondria (matrix) | function: mitochondrial import targeting (cleaved)
- binding partner: TOM20 import receptor | surface: TOM20 hydrophobic groove (phi-chi-chi-phi-phi)
- AKAP-like because: amphipathic with hydrophobic + charged face
- NOT AKAP because: Arg-rich, acidic-depleted targeting peptide; cleaved; recognized by TOM/TIM not PKA
- context evidence: classic cleavable MTS; structural: TOM20-presequence complexes (Saitoh 2007)
- expected failure mode: likely PSSM<12 (different anchors) -> PSSM filters
- intended use: **stress_test** | priority: **medium** | verification: verify_region
- notes: if PSSM<12, report PSSM already filters this class

### MT02 — Ornithine carbamoyltransferase presequence (Homo sapiens)
- accession: **P00480** (UniProt); PDB: —; region: N-terminal presequence ~1-32; type: full_protein
- location: mitochondria (matrix) | function: import targeting (cleaved)
- binding partner: TOM20 | surface: TOM20 groove
- AKAP-like because: amphipathic targeting helix
- NOT AKAP because: positive-inside, acidic-poor; TOM/TIM substrate
- context evidence: canonical OTC presequence; structural: —
- expected failure mode: likely PSSM<12
- intended use: **stress_test** | priority: **low** | verification: verify_region

### MT03 — Cytochrome c oxidase subunit 4 isoform 1 presequence (Homo sapiens)
- accession: **P13073** (UniProt); PDB: —; region: N-terminal presequence ~1-22; type: full_protein
- location: mitochondria | function: import targeting (cleaved)
- binding partner: TOM20 | surface: TOM20 groove
- AKAP-like because: amphipathic presequence helix
- NOT AKAP because: Arg-rich acidic-poor targeting peptide
- context evidence: model MTS (COX4); structural: —
- expected failure mode: likely PSSM<12
- intended use: **stress_test** | priority: **low** | verification: verify_region

## Class 4 — Lipid-binding / apolipoprotein helices

### LP01 — Apolipoprotein A-I (Homo sapiens)
- accession: **P02647** (UniProt); PDB: 1AV1; region: class-A 11/22-mer repeats (mature); type: full_protein
- location: secreted (HDL) | function: lipid binding / HDL scaffold
- binding partner: phospholipid / lipoprotein surface | surface: lipid-facing hydrophobic face
- AKAP-like because: textbook class-A amphipathic helix, sharp face separation, high moment
- NOT AKAP because: lipid binder; interfacial-basic + center-acidic snorkel; tandem repeats; not PKA
- context evidence: major HDL apolipoprotein; structural: PDB 1AV1 (lipid-free belt)
- expected failure mode: HIGH danger-zone entry; central acidic cluster may trigger negdet->unlikely
- intended use: **training_candidate+stress_test** | priority: **high** | verification: verify_region
- notes: also tests whether negdet rule earns its keep

### LP02 — Apolipoprotein E (C-terminal lipid-binding) (Homo sapiens)
- accession: **P02649** (UniProt); PDB: 1B68; region: C-terminal lipid-binding helices ~244-272; type: full_protein
- location: secreted | function: lipoprotein lipid binding
- binding partner: lipoprotein lipid surface | surface: lipid-facing face
- AKAP-like because: amphipathic lipid-binding helices
- NOT AKAP because: lipid target; not PKA D/D
- context evidence: ApoE lipid binding; structural: PDB 1B68 region
- expected failure mode: danger-zone entry plausible
- intended use: **stress_test** | priority: **medium** | verification: verify_region

### LP03 — Perilipin-1 (class-A repeats) (Homo sapiens)
- accession: **O60240** (UniProt); PDB: —; region: 11-mer amphipathic repeat region; type: full_protein
- location: lipid droplet | function: lipid droplet coat
- binding partner: lipid droplet monolayer | surface: lipid-facing face
- AKAP-like because: class-A amphipathic 11-mer repeats
- NOT AKAP because: lipid-droplet binder; not PKA (note: PLIN1 is regulated by PKA but anchoring helix is lipid-facing)
- context evidence: PAT-family LD protein; structural: —
- expected failure mode: danger-zone entry plausible
- intended use: **stress_test** | priority: **medium** | verification: verify_region
- notes: PLIN1 is a PKA SUBSTRATE - ensure region chosen is the lipid-binding helix, not a PKA site

## Class 5 — Generic helix-in-groove PPI motifs

### PG01 — BAD BH3 helix (Homo sapiens)
- accession: **Q92934** (UniProt); PDB: 1G5J; region: BH3 ~103-127; type: full_protein
- location: cytosol/mitochondria | function: apoptosis regulation (pro-apoptotic BH3-only)
- binding partner: Bcl-XL / Bcl-2 hydrophobic groove | surface: protein-groove interface
- AKAP-like because: short amphipathic helix docking into a hydrophobic protein groove (mechanistically closest)
- NOT AKAP because: BH3 anchor spacing (Lxxxx(D/phi)) and Bcl-2 partner differ from AKAP D/D
- context evidence: BH3-in-groove binding; structural: PDB 1G5J
- expected failure mode: anchor spacing differs; possible PSSM<12
- intended use: **stress_test** | priority: **high** | verification: verify_region
- notes: most diagnostic if it reaches danger zone

### PG02 — BIM BH3 helix (Homo sapiens)
- accession: **O43521** (UniProt); PDB: 1PQ1; region: BH3 ~141-166; type: full_protein
- location: cytosol/mitochondria | function: apoptosis (BH3-only)
- binding partner: Bcl-2 family groove | surface: protein-groove interface
- AKAP-like because: helix-in-groove amphipathic
- NOT AKAP because: Bcl-2 partner; BH3 motif spacing
- context evidence: BH3 binding; structural: PDB 1PQ1
- expected failure mode: anchor spacing differs
- intended use: **stress_test** | priority: **medium** | verification: verify_region

### PG03 — p53 transactivation helix (Homo sapiens)
- accession: **P04637** (UniProt); PDB: 1YCR; region: TAD helix ~17-29; type: full_protein
- location: nucleus | function: transcription; MDM2 regulation
- binding partner: MDM2 hydrophobic cleft | surface: protein-groove interface
- AKAP-like because: short helix in hydrophobic groove (F19/W23/L26 anchors)
- NOT AKAP because: MDM2 partner; 3-residue anchor spacing; very short (~3 turns)
- context evidence: p53-MDM2 (Kussie 1996); structural: PDB 1YCR
- expected failure mode: too short for 5-turn AKAP support; likely filtered
- intended use: **stress_test** | priority: **high** | verification: verify_region
- notes: tests 5-turn support discriminator

### PG04 — NCOA1/SRC-1 LxxLL NR box (Homo sapiens)
- accession: **Q15788** (UniProt); PDB: 1GWQ; region: NR box LxxLL motif; type: full_protein
- location: nucleus | function: nuclear receptor coactivation
- binding partner: nuclear receptor LBD AF-2 groove | surface: protein-groove interface
- AKAP-like because: short amphipathic helix in a hydrophobic groove
- NOT AKAP because: LxxLL motif; NR-LBD partner; short
- context evidence: NR box binding; structural: PDB 1GWQ
- expected failure mode: short motif; likely filtered
- intended use: **stress_test** | priority: **medium** | verification: verify_region

## Class 6 — DDIP-like / DPY30-binding helices

### DD01 — DPY30 (Homo sapiens)
- accession: **Q9C005** (UniProt); PDB: 3G36; region: dimerization/binding region; type: full_protein
- location: nucleus | function: MLL/COMPASS assembly; binds D/D-like fold
- binding partner: ASH2L (Sdc-binding) / D/D-like domain | surface: D/D-like domain interface (NOT PKA R-subunit)
- AKAP-like because: engages a dimerization/docking (D/D-like) fold via helix
- NOT AKAP because: partner is DPY30/ASH2L system, not PKA RI/RII; architecture differs
- context evidence: COMPASS subunit; structural: PDB 3G36 region
- expected failure mode: D/D-like but non-PKA; tests DDIP discrimination
- intended use: **stress_test** | priority: **high** | verification: verify_region_and_AKAP_status
- notes: FLAGGED: confirm it is experimentally NOT a PKA anchor before training use

### DD02 — ASH2L (DPY30-binding helix) (Homo sapiens)
- accession: **Q9UBL3** (UniProt); PDB: —; region: DPY30-binding region; type: full_protein
- location: nucleus | function: COMPASS H3K4 methyltransferase complex
- binding partner: DPY30 D/D-like domain | surface: D/D-like interface
- AKAP-like because: helix binding a D/D-like fold
- NOT AKAP because: non-PKA D/D-like partner
- context evidence: ASH2L-DPY30 interaction; structural: —
- expected failure mode: DDIP-like discrimination test
- intended use: **stress_test** | priority: **medium** | verification: verify_region_and_AKAP_status
- notes: FLAGGED ambiguous: confirm non-PKA-anchoring before any training use

## Class 7 — PDE / GAF-domain regulatory & dimerization amphipathic helices (parent: generic_intramolecular_or_dimerization_amphipathic_helix)

### PDE01 — PDE2A (GAF-B regulatory helix, 404-427) (Homo sapiens)
- accession: **O00408** (UniProt); PDB: 3IBJ; region: 404-427 (core VSVLLQEIITEA); GAF-B domain ~393-541; type: full_protein
- location: cytosol / mitochondria-associated (isoform-dependent) | function: cyclic-nucleotide phosphodiesterase regulatory (GAF) domain; cGMP-activated
- binding partner: cGMP (GAF-B allosteric site) / PDE2A homodimer interface / intramolecular domain packing | surface: regulatory-domain helix; domain-packing / dimerization-associated face (NOT PKA D/D groove)
- AKAP-like because: amphipathic helix, high RII PSSM (15.14), ML v2 high score (0.992)
- NOT AKAP because: PDE2A is a phosphodiesterase; window lies in the GAF-B regulatory domain, not a known PKA RI/RII D/D anchoring domain; hydrophobic face supports domain packing / dimerization, not PKA anchoring
- context evidence: PDE2A domains N-term 1-214, GAF-A 215-372, GAF-B 393-541, catalytic 579-941; homodimer (PNAS 2009, 0907635106); structural: PDB 3IBJ: GAF/catalytic domains joined by long alpha-helices; dimer interface spans the molecule
- expected failure mode: high-PSSM amphipathic non-AKAP regulatory/domain-packing helix; promoted via high-bg downgrade (pssm>=14 AND ml>=0.95)
- intended use: **training_candidate+stress_test** | priority: **very_high** | verification: verify_accession_region
- notes: CONFIRMED AKAPSpred v5.1 + ML v2 danger-zone FALSE POSITIVE (final tier=high). CLASSIFICATION: biologically supported non-AKAP contextual false positive (domain/function annotation); NOT experimentally proven non-binder — no direct PKA RI/RII binding assay exists. Burden of proof is on the AKAP-positive call.

### PDE02 — PDE3A (N-terminal region, 64-87) (Homo sapiens)
- accession: **Q14432** (UniProt); PDB: —; region: 64-87 (core LSFLLALLVRLV); type: full_protein
- location: membrane-associated / cytosol | function: cGMP-inhibited cyclic-nucleotide phosphodiesterase
- binding partner: lipid membrane / intramolecular (N-terminal hydrophobic region) | surface: membrane-association / domain context (NOT PKA D/D groove)
- AKAP-like because: hydrophobic N-terminal region scored by sensitive PSSM scan
- NOT AKAP because: PDE3A is a phosphodiesterase; window not a known PKA D/D anchoring helix
- context evidence: UniProt Q14432 phosphodiesterase; structural: —
- expected failure mode: below danger zone (pssm 7.64<12) — should be filtered
- intended use: **same_family_negative_control** | priority: **medium** | verification: verify_region
- notes: CORRECTLY REJECTED by v5.1+ML v2 (sensitive_only; pssm 7.64, ml 0.005, amphipathic=False). All three layers reject. Same-family negative control.

### PDE03 — PDE4D (267-290) (Homo sapiens)
- accession: **Q08499** (UniProt); PDB: —; region: 267-290 (core YQKLASETLEEL); type: full_protein
- location: cytosol / membrane (isoform-dependent) | function: cAMP-specific cyclic-nucleotide phosphodiesterase
- binding partner: recruited to AKAP/signalosome complexes (e.g. mAKAP, AKAP9) but not itself a PKA anchor | surface: regulatory/UCR domain context (NOT PKA D/D groove)
- AKAP-like because: amphipathic window scored by sensitive PSSM scan
- NOT AKAP because: PDE4D is a phosphodiesterase recruited BY AKAPs, not a PKA-anchoring AKAP itself
- context evidence: UniProt Q08499; PDE4D is a signalosome effector; structural: —
- expected failure mode: borderline; rejected (pssm 11.04<12, ml 0.272)
- intended use: **same_family_negative_control** | priority: **medium** | verification: verify_accession_region
- notes: CORRECTLY REJECTED by v5.1+ML v2 (sensitive_only; pssm 11.04, ml 0.272). Same-family negative control. NB: PDE4D is recruited by AKAPs but is not an AKAP.

### PDE04 — PDE2A (302-325, negative-determinant control) (Homo sapiens)
- accession: **O00408** (UniProt); PDB: —; region: 302-325 (core LKDLTSEDVQQL); GAF-A domain ~215-372; type: full_protein
- location: cytosol / mitochondria-associated | function: cyclic-nucleotide phosphodiesterase regulatory (GAF-A) domain
- binding partner: GAF-A dimerization locus / intramolecular packing | surface: regulatory-domain helix (NOT PKA D/D groove)
- AKAP-like because: amphipathic window in a GAF regulatory domain
- NOT AKAP because: phosphodiesterase GAF-A region; also carries a negative determinant on the face
- context evidence: PDE2A GAF-A 215-372 (dimerization locus, PNAS 2002); structural: —
- expected failure mode: rejected by negative-determinant red flag (n_negdet=1)
- intended use: **same_protein_negative_determinant_control** | priority: **medium** | verification: verify_region
- notes: CORRECTLY REJECTED by v5.1+ML v2 (unlikely; biological red flag, n_negdet=1). Same-protein control showing the negative-determinant path firing.

