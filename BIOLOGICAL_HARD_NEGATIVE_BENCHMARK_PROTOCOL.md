# Biological Hard-Negative Benchmark Protocol — AKAPSpred v5.1 + ML v2

**Purpose.** Test whether the deployed pipeline distinguishes true PKA-anchoring
AKAP helices from *other real amphipathic helices* that have different cellular
functions, binding partners, local context, and structural constraints. This is
a diagnostic, not a training run. No retraining, no threshold change, no model
replacement.

## Central biological idea

A true AKAP anchoring helix is not merely amphipathic. It binds the PKA RI/RII
D/D domain, presents a hydrophobic face compatible with that groove, spans enough
turns (often ~5), avoids acidic negative determinants on the binding face, and sits
in a protein context compatible with scaffolding/localization. Many non-AKAP
helices are amphipathic, but their hydrophobic face docks into something else — a
partner helix, a lipid bilayer, a mitochondrial import receptor, or a different
protein groove. Each hard negative is therefore selected and interpreted by its
**real biological role**, recorded in mandatory context fields, not by sequence
pattern alone.

## The discriminating principle (what separates AKAP from each class)

The shared feature is the amphipathic-helix *form*; the differences are (a) what
the hydrophobic face engages and (b) the composition/charge of the polar face. The
polar face is the most class-specific signal and is the leading candidate
discriminator the current model does not explicitly encode.

## Biological hard-negative classes

Full per-candidate detail (accession, PDB, region, partner, surface, evidence,
expected failure mode) is in `hard_negative_candidate_list_biological_context.md`
and `.csv`. Summary:

1. **Coiled-coil** — hydrophobic seam packs against a partner helix (heptad
   3.5 res/turn); discriminator: coiled-coil/heptad periodicity, long repeat
   continuity, helix–helix packing context. *training_candidate + stress_test.*
2. **Membrane-binding amphipathic** (ALPS, α-synuclein, epsin H0, AMPs) —
   hydrophobic face inserts into lipid; polar face is charge-poor (ALPS) or
   cationic (AMPs); discriminator: polar-face charge pattern, interfacial
   aromatics. *Highest-priority class; training_candidate + stress_test.*
3. **Mitochondrial targeting presequences** — Arg-rich, acidic-depleted,
   TOM20-recognized, cleaved; discriminator: N-terminal position, high Arg/low
   D/E, TOM20 motif. Often filtered by PSSM. *stress_test first.*
4. **Lipid-binding / apolipoprotein** (ApoA-I, ApoE, perilipin) — class-A helix
   with interfacial-basic / center-acidic snorkel, tandem 11/22-mers; the central
   acidic cluster may trip the AKAP negative-determinant rule. *training_candidate
   + stress_test.*
5. **Generic helix-in-groove PPI** (BH3, p53–MDM2, LxxLL) — mechanistically
   closest (helix in a hydrophobic groove) but partner and anchor spacing differ
   and many are too short for 5-turn support. *high-value stress_test.*
6. **DDIP-like / DPY30-binding** — engage a D/D-*like* fold but not PKA RI/RII;
   discriminator: shorter helix, non-PKA partner specificity. *stress_test; only
   train if non-PKA-anchoring status is experimentally confirmed.*

## Controls — kept strictly separate (via `set_type`)

- **positive_control** — known RI/RII/dual AKAP anchoring helices; confirm the
  pipeline still detects true AKAPs.
- **easy_negative_control** — random / composition-shuffled; sanity only.
- **synthetic_mechanistic_decoy** — mutated AKAP motifs (anchor disruption,
  Pro/Gly breaker, D/E on the hydrophobic face, shortened helix); mechanistic
  stress only.

Synthetic decoys must **not** dominate or substitute for real biological hard
negatives. QC enforces that controls are not mixed with biological hard negatives.

## Danger-zone logic

Primary diagnostic subset:

```
label = non_AKAP   AND   amphipathic = True   AND   pssm_score >= 12
```

For every danger-zone hit, the harness records and the report interprets:

1. full protein or peptide-only control? (`is_peptide_only`)
2. is the predicted window in the known functional helix region? (`in_known_helix_region`)
3. is the known partner PKA RI/RII or something else? (`partner_is_PKA` = no by construction; `known_binding_partner` shown)
4. what does the hydrophobic face bind — protein/lipid/membrane/another helix? (`known_binding_surface`)
5. is the local cellular context compatible with AKAP function? (`cellular_location`)
6. real biological near-miss or artificial decoy? (`set_type`)

## Interpretation rules

- Class with **no** danger-zone entry → PSSM already filters this biological class.
- Danger-zone entry but **low** ML probability → ML v2 gives useful independent
  discrimination for that context.
- Danger-zone entry **and high** ML probability → priority hard-negative class for
  ML v3.
- FPs mainly from **peptide-only** examples → report separately (no protein
  background context).
- FPs from **full proteins in their known helix region** → strongest evidence of a
  real deployment-relevant weakness.

## ML v3 implications (per class)

| Class | v3 role | Recommended v3 feature(s) |
|---|---|---|
| Coiled-coil | training_candidate | coiled-coil/heptad periodicity; helix length / 5-turn-plus continuity; hydrophobic-face continuity |
| Membrane amphipathic | training_candidate | polar-face charge distribution (cationic / charge-poor); interfacial aromatic enrichment; hydrophobic-face continuity |
| Mitochondrial targeting | stress_test_only (train only if danger-zone breaches) | N-terminal targeting-presequence pattern; high Arg / low D/E |
| Apolipoprotein / lipid | training_candidate | lipid-binding snorkel pattern; central polar-face acidic cluster; tandem 11/22-mer periodicity; hydrophobic-face charge disruption |
| Helix-in-groove PPI | stress_test_only (holdout) | helix-in-groove anchor spacing; 5-turn AKAP support; helix length |
| DDIP-like | stress_test_only until labels clean | DDIP-like short-helix support; D/D-like-but-non-PKA partner signature |

Cross-cutting v3 features motivated by this benchmark: hydrophobic-face
continuity, 5-turn AKAP support, hydrophobic-face charge disruption, polar-face
charge distribution. Easy negatives are already handled by PSSM, so v3's real
negatives should concentrate in the **amphipathic + high-PSSM** region.

## Sourcing verified sequences (no fabrication)

Sequences are added only after verification. Two routes:

1. **Curator pastes/fetches** each accession's sequence (UniProt
   `https://rest.uniprot.org/uniprotkb/{ACC}.fasta`, or PDB SEQRES) and records it
   with `fetched_from_url`, `fetch_date`, `verified_by` in the metadata.
2. **Fetch by accession** from the candidate list once each `verification_status`
   is cleared (every row currently says `verify_*`; confirm accession, isoform,
   and region first — several regions are approximate).

FASTA header convention: `>{id}` (e.g. `>MB02`) matching the metadata `id`, so QC
and the harness join cleanly. Full proteins are preferred (realistic
`background_risk`); peptide-only entries are allowed but analysed separately.

## Workflow

```
# 1. clear verification_status, fill sequences -> hard_negative_amphipathic_set.fasta
#    and a verified metadata CSV (use the template / candidate list)
python3 hard_negative_qc.py \
    --meta  hard_negative_metadata_biological_context.csv \
    --fasta hard_negative_amphipathic_set.fasta          # must PASS

python3 run_hard_negative_benchmark.py \
    --fasta hard_negative_amphipathic_set.fasta \
    --meta  hard_negative_metadata_biological_context.csv \
    --screen ../akap_screen.py --outdir .
# -> hard_negative_v51_screening.csv
#    hard_negative_danger_zone.csv
#    hard_negative_class_summary.csv
#    HARD_NEGATIVE_DIAGNOSTIC_REPORT.md
```

## Constraints honoured
No fabricated sequences. No ML v3 retraining. No v5.1 threshold change. No ML v2
replacement. Synthetic decoys are mechanistic-only and segregated. Ambiguous-AKAP
candidates are flagged (QC blocks their use as training negatives) rather than
silently treated as negatives.
