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

## Beginner manual

For non-expert users, see **`USER_MANUAL.md`**. It explains:

- Which input mode to use.
- How to fetch proteins from UniProt.
- What each result column means.
- How to decide whether a hit is strong or weak.
- What to do after getting predicted hits.

### Fast interpretation guide

| Result item | Good sign | Caution sign |
|---|---|---|
| `classification` | `AKAP` | `ambiguous`, `DDIP`, or `unlikely` |
| `pssm_score` | Higher than other hits of the same isoform | Near threshold |
| `ml_prob` | Close to 1.0 | Around or below 0.5 |
| `n_negdet` | 0 | 1 or more |
| `amphipathic` | True | False |
| `canonical` | True supports the hit | False is not automatic rejection |

Use the screener as a prioritization tool. A predicted motif is a candidate for follow-up, not final proof of PKA binding.

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
