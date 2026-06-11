# 🧬 AKAP Domain Screener

Screen proteins for **A-Kinase Anchoring Protein (AKAP)** amphipathic-helix motifs that bind PKA regulatory subunits.

Based on: **Burgers et al. (2015)** *"A Systematic Evaluation of Protein Kinase A–A-Kinase Anchoring Protein Interaction Motifs"*, Biochemistry 54, 11–21. [DOI: 10.1021/bi500721a](https://doi.org/10.1021/bi500721a)

---

## What it does

An AKAP domain is a short amphipathic α-helix (3–4 turns) that docks its hydrophobic face onto the D/D domain of the PKA regulatory-subunit dimer. This tool screens protein sequences for these motifs using two complementary methods:

| Method | Engine | Performance |
|--------|--------|-------------|
| **PSSM profile scan** | Log-odds matrix from 849 RIIα + 28 RIα motifs | RII: 99.2% recall, 0.8% FPR · RI: 100% recall |
| **ML classifier** | GradientBoosting on biophysical features | RII: AUC 0.999 · RI: AUC 1.000 |

**Features detected by the ML model** (ranked by importance):
1. PSSM positional score
2. Helix propensity (Pace-Scholtz)
3. Hydrophobic contrast (anchor vs polar face)
4. Minimum anchor hydrophobicity
5. FFT amphipathicity at α-helix period (3.6 residues)

---


---

## Beginner guide for non-expert users

This section explains how to use the app and how to interpret the values in plain language. The screener is a prioritization tool: a predicted motif is a candidate for follow-up, not final proof of PKA binding.

## 1. What this tool does

The AKAP Domain Screener searches protein sequences for short regions that look like **AKAP motifs**.

An **AKAP motif** is usually a short amphipathic alpha helix. In simple terms, this means one side of the helix is hydrophobic and can dock onto the dimerization/docking domain of PKA regulatory subunits.

The tool is useful for:

- Finding possible AKAP-like motifs in a protein.
- Comparing candidate motifs across many proteins.
- Prioritizing motifs for experimental validation.
- Deciding which short peptide regions may be worth testing.

The tool is **not** final proof of PKA binding. It gives candidate regions that need biological interpretation and validation.

---

## 2. How to run the app

Open a terminal in the folder containing the files and run:

```bash
pip install -r requirements.txt
streamlit run akap_app.py
```

A browser window should open automatically.

---

## 3. Easiest use mode: search UniProt

For most users, use **Search by protein name / UniProt ID**.

### Steps

1. Select **Search by protein name / UniProt ID**.
2. Type one protein, gene name, or UniProt accession per line.

Example:

```text
Ezrin
AKAP10
WAVE1
P15311
O43572
```

3. Choose organism, usually **Human**.
4. Keep **Reviewed (Swiss-Prot) only** checked for cleaner results.
5. Click **Search UniProt**.
6. Check the fetched protein table.
7. Click **Run AKAP Screen**.

### Important note

If you type a common protein name, UniProt may return several proteins. The app automatically uses the **top hit**. Before using results for a manuscript, check that the accession, organism, and sequence length are correct.

---

## 4. Other input modes

### Paste sequence(s)

Use this when you already have the protein sequence.

FASTA example:

```text
>MyProtein
MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDW...
```

Raw sequence is also allowed, but FASTA is better because it keeps the protein name.

### Upload FASTA file

Use this for one or many protein sequences saved as `.fasta`, `.fa`, `.faa`, or `.txt`.

### Upload CSV file

Use this for a table of proteins.

Accepted column names are flexible:

| Purpose | Column names that work |
|---|---|
| Protein name | `protein`, `name`, `gene`, `label` |
| UniProt ID | `uniprot`, `accession`, `entry` |
| Sequence | `sequence`, `seq` |

CSV example:

```csv
protein_name,uniprot_id,sequence
AKAP10,O43572,
Ezrin,,
My_custom_protein,,MAAADSGRLHAAAL...
```

If a sequence is missing, the app tries to fetch it from UniProt.

---

## 5. Recommended settings for non-experts

Use the default settings first.

| Setting | Recommended value | Meaning |
|---|---:|---|
| PKA-RIIα threshold | 7.0 | Sensitive default for RII-like motifs |
| PKA-RIα threshold | 12.0 | Sensitive default for RI-like motifs |
| Require literal regex match | OFF | Leave off; turning it on may miss real motifs |
| Enable ML scoring | ON if available | Adds an extra machine-learning confidence score |
| ML probability threshold | 0.5 | Balanced default |

For stricter screening, increase thresholds or set ML probability to 0.7–0.8.

---

## 6. What each result value means

### Main columns

| Column | Simple meaning | How to interpret |
|---|---|---|
| `protein` | Protein name or accession | Which protein contains the predicted motif |
| `isoform` | RII or RI | Which PKA regulatory subunit type the motif resembles |
| `classification` | Overall motif class | Start here first |
| `start`, `end` | Position in the protein | Location of the predicted motif |
| `core` | Main motif sequence | Short region most important for docking |
| `window` | Longer scored sequence | Full sequence window used by the model |
| `pssm_score` | Similarity to known motifs | Higher is stronger; rank RI and RII separately |
| `ml_prob` | ML confidence | Near 1.0 is stronger; around 0.5 is moderate |
| `dual` | Overlapping RI/RII prediction | May indicate a region compatible with both RI and RII |
| `canonical` | Strict motif-pattern match | Helpful if true, but false does not automatically mean wrong |
| `amphipathic` | Helix face pattern check | True supports AKAP-like behavior |
| `pI` | Estimated charge property | Annotation only; not pass/fail |
| `n_negdet` | Number of negative determinants | 0 is best; more can weaken confidence |
| `negdet_severity` | Severity of negative determinants | `none` is best; higher severity means more caution |
| `helix_turns` | Approximate hydrophobic-face length | More turns generally support AKAP-like binding |
| `dd_class` | Predicted D/D-domain interaction class | Helps separate AKAP-like vs DDIP-like behavior |

---

## 7. How to interpret classification

| Classification | Meaning | Priority |
|---|---|---|
| **AKAP** | Strong candidate PKA-anchoring motif | High priority |
| **ambiguous** | Possible motif but not fully clear | Medium priority; inspect manually |
| **DDIP** | May bind a D/D-like groove but not predicted as strong PKA anchor | Context-dependent |
| **unlikely** | Weak or disrupted motif | Low priority |

---

## 8. What is a good hit?

A stronger candidate usually has:

- `classification = AKAP`
- High `pssm_score`
- High `ml_prob`, if ML is available
- `n_negdet = 0`
- `amphipathic = True`
- A reasonable helix-like hydrophobic pattern in the helical wheel plot
- Biological context that makes sense for PKA regulation

A weaker candidate may have:

- `classification = ambiguous` or `unlikely`
- Low PSSM score
- Low ML probability
- Many negative determinants
- Acidic residues on hydrophobic anchor positions

---

## 9. What to do after getting hits

For publication-quality analysis, do not stop at the table. Recommended follow-up:

1. Confirm the protein accession and isoform.
2. Check whether the predicted region is conserved.
3. Check secondary-structure tendency with AlphaFold, NetSurfP, PSIPRED, or similar tools.
4. Compare with known AKAP motifs in the literature.
5. Design mutations at hydrophobic anchor positions.
6. Test binding experimentally, for example peptide binding, pull-down, co-IP, or localization assay.

---

## 10. Common problems

### “I clicked Search UniProt but Run AKAP Screen shows nothing.”

Use the updated `akap_app.py`. It saves fetched UniProt proteins in Streamlit session state so they remain available after the next button click.

### “No hit was found.”

This means no region passed the current thresholds. It does not prove the protein cannot interact with PKA. Check the sequence, use default thresholds, and consider splice isoforms.

### “Too many hits were found.”

Increase thresholds or set ML probability threshold to 0.7–0.8. Then prioritize `classification = AKAP` and `n_negdet = 0`.

### “The protein name gives the wrong UniProt entry.”

Use the exact UniProt accession instead of a gene/protein name.

---

## 11. One-sentence interpretation template

Use this wording when reporting a result:

> The AKAP Domain Screener predicted an AKAP-like motif in [protein] at residues [start]–[end], with [isoform] preference, PSSM score [score], ML probability [probability], and [number] negative determinants. This region should be considered a candidate motif for experimental validation.

---

## Quick start

### Web app (Streamlit)

```bash
pip install -r requirements.txt
streamlit run akap_app.py
```

The app accepts:
- **Paste** a sequence or FASTA
- **Upload** a FASTA or CSV file
- **Search by protein name** — type `Ezrin`, `AKAP10`, `WAVE1` etc. and the app fetches sequences from UniProt automatically

### Command line

```bash
# Screen a FASTA file
python akap_screen.py proteins.fasta -o hits.csv

# Screen from a CSV (protein names / UniProt IDs / sequences)
python akap_from_csv.py my_list.csv -o hits.csv

# ML predictions on a FASTA
python akap_ml.py --predict proteins.fasta -o ml_hits.csv
```

---

## Installation

```bash
git clone https://github.com/<your-username>/akap-screener.git
cd akap-screener
pip install -r requirements.txt
```

### Retrain the ML model (optional)

The pre-trained model (`akap_ml_model.joblib`) is included. To retrain from the SI data:

```bash
# Place the SI xlsx files in the same directory, then:
python akap_ml.py
```

---

## Files

| File | Description |
|------|-------------|
| `akap_app.py` | 🌐 Streamlit web app — main entry point for interactive use |
| `akap_screen.py` | 🔬 Core PSSM screening engine (CLI) |
| `akap_pssm.json` | 📊 Position-specific scoring matrix (849 RII + 28 RI motifs) |
| `akap_ml.py` | 🤖 ML classifier — train, evaluate, predict |
| `akap_ml_model.joblib` | 💾 Pre-trained GradientBoosting model |
| `akap_from_csv.py` | 📋 CLI wrapper: CSV → UniProt fetch → screen |
| `validate_akap_screen.py` | ✅ Validation script against SI data |
| `requirements.txt` | 📦 Python dependencies |

### Supporting data (not included — obtain from the publisher)

| File | Source |
|------|--------|
| `bi500721a_si_001.xlsx` | SI Table 1: 849 PKA-RIIα motifs |
| `bi500721a_si_002.xlsx` | SI Table 2: 28 PKA-RIα motifs |

These are needed only for retraining (`akap_ml.py`) and validation (`validate_akap_screen.py`). The pre-built PSSM and ML model already encode this data.

---

## Input formats

### FASTA

```
>AKAP10
MAAADSGRLH...EAQEELAWKIAKMIVSDIMQQAQY...
>MyProtein
MKWVTFIS...
```

### CSV

Any CSV/TSV with columns matching these keywords (case-insensitive, order doesn't matter):

| Column keywords | What it's used for |
|---|---|
| `protein`, `name`, `gene` | Display name |
| `uniprot`, `accession` | UniProt ID → auto-fetch sequence |
| `sequence`, `seq` | Amino acid sequence (used directly) |

```csv
protein_name,uniprot_id,sequence
AKAP10,O43572,
Ezrin,,
My_custom_protein,,MAAADSGRLH...
```

The tool resolves missing sequences by: sequence column → UniProt ID fetch → protein name search.

---

## Output

| Column | Description |
|--------|-------------|
| `protein` | Protein identifier |
| `isoform` | `RII` or `RI` — which PKA-R isoform |
| `dual` | `True` if both RI and RII sites overlap → dual-specific AKAP |
| `win_start`, `win_end` | Window position in the protein |
| `core` | The 12-mer (RII) or 18-mer (RI) hydrophobic core |
| `window` | Full 24-mer (RII) or 30-mer (RI) aligned window |
| `pssm_score` | Log-odds score — higher = more AKAP-like |
| `ml_prob` | ML classifier probability (0–1) |
| `canonical` | Does the core match the strict consensus regex? |
| `amphipathic` | Passes the hydrophilic-face polarity check? |
| `pI` | Isoelectric point of the core (annotation) |
| `helix_approx` | Approximate helix propensity (0–1) |
| `classification` | **AKAP** / **DDIP** / **ambiguous** / **unlikely** (Falcone & Scott 2025) |
| `n_negdet` | Count of charged residues (D/E) at hydrophobic anchor positions |
| `negdet_severity` | `none` / `mild` (1) / `severe` (2+) |
| `contiguous_hydro_turns` | Estimated amphipathic helix turns (AKAP ≥5, DDIP 3–4) |
| `dd_class` | Predicted d/d domain partner: `RIID2` / `RID2` / `DPY-30` |

---

## How the PSSM thresholds were calibrated

| Isoform | Default threshold | Recall on SI data | Null FP rate |
|---------|-------------------|-------------------|--------------|
| PKA-RIIα | 7.0 | 99.2% (842/849) | 0.8% |
| PKA-RIα | 12.0 | 100% (28/28) | 0.4% |

The 7 missed RIIα motifs are the paper's own weakest hits (high MAST E-values). Lowering to 6.0 gives 100% recall at ~1.1% FPR.

---

## Citation

If you use this tool, please cite the original papers:

> Burgers PP, van der Heyden MAG, Kok B, Heck AJR, Scholten A. (2015)
> A Systematic Evaluation of Protein Kinase A–A-Kinase Anchoring Protein Interaction Motifs.
> *Biochemistry* 54(1):11–21. [doi:10.1021/bi500721a](https://doi.org/10.1021/bi500721a)

> Falcone JI, Scott JD. (2025)
> The ascent of AKAPs, from architectural elements to kinase anchors: a perspective.
> *Biochem J* 482(10):485–498. [doi:10.1042/BCJ20253085](https://doi.org/10.1042/BCJ20253085)

---

## DDIP vs AKAP classification (NEW)

Based on Falcone & Scott (2025), the tool now distinguishes true PKA-anchoring AKAPs from
D/D domain interacting proteins (DDIPs) — proteins that bind the same d/d groove via amphipathic
helices but do NOT anchor PKA.

| Classification | Meaning | Criteria |
|---|---|---|
| **AKAP** | Strong PKA anchor | ≥5 helix turns, no negative determinants, high PSSM |
| **DDIP** | D/D domain interactor (not PKA) | 3–4 helix turns (shorter amphipathic helix) |
| **ambiguous** | Could be either | Borderline; needs experimental validation |
| **unlikely** | Disrupted hydrophobic face | Charged residue (D/E) at anchor position |

### Negative determinants

A single Asp or Glu on the hydrophobic face of the amphipathic helix abolishes PKA binding.
This was demonstrated in:
- **OPA1 fungal forms** (Asp at position 6 → no AKAP function)
- **smAKAP S66D mutant** (Ser→Asp at anchor → loss of PKA binding)

The tool flags these as `negdet_severity: mild (1 violation) or severe (2+)`.

### d/d domain class

Each hit is annotated with its predicted binding partner class:
- **RIID2** — PKA-RIIα-like d/d domain (most AKAPs)
- **RID2** — PKA-RIα-like d/d domain (dual-specific AKAPs)
- **DPY-30** — histone methylation / dosage compensation machinery

---

## License

MIT


---

🧬 AKAP Domain Screener — kowith@ccs.tsukuba.ac.jp


## v4 update: DDIP/AKAP and RI/RII specificity

This version adds a practical interpretation layer inspired by Falcone & Scott (2025):

- **classification**: AKAP, DDIP, ambiguous, or unlikely.
- **n_negdet / negdet_severity**: flags Asp/Glu residues on hydrophobic anchor positions. These are negative determinants for PKA anchoring.
- **helix_turns**: estimated amphipathic-helix extent. Longer helices support AKAP-like anchoring; shorter 3–4 turn helices may be DDIP-like.
- **dd_class**: predicted d/d-domain partner class, currently RIID2 for RII-like hits and RID2 for RI-like hits.
- **predicted_specificity**: full-protein RI/RII preference call. This is important for proteins such as **SPHKAP**, which should be interpreted as **RI-specific** even if a weaker RII-like window also crosses the RII threshold.
- **ri_best / rii_best / ri_rii_ratio**: the best RI score, best RII score, and their ratio across the full sequence.

For non-expert users: prioritize hits with **classification = AKAP**, **n_negdet = 0**, and a specificity label that matches your biological question. Treat **DDIP**, **ambiguous**, and **unlikely** as lower-confidence candidates that need experimental validation.
