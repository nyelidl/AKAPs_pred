#!/usr/bin/env python3
"""
akap_app.py — Streamlit web app for AKAP domain screening.

Run:  streamlit run akap_app.py
"""

import io
import json
import math
import os
import re
import sys
import time
from collections import Counter

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# ─── Ensure our modules are importable ───
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

import akap_screen as A

# ─── Load ML model bundle (prefer v2) ───
_HAVE_ML = False
_ML_BUNDLE = None
_ML_MODEL_PATH = None
_ML_METADATA = {}
try:
    # Use backend loader so CLI and app share the same model loading behavior.
    _ML_BUNDLE = A.load_ml_model()
    _HAVE_ML = _ML_BUNDLE is not None
    for _candidate in ("akap_ml_model_v2.joblib", "akap_ml_model.joblib"):
        _p = os.path.join(SCRIPT_DIR, _candidate)
        if os.path.exists(_p):
            _ML_MODEL_PATH = _p
            break
    _meta_path = os.path.join(SCRIPT_DIR, "akap_ml_model_v2_metadata.json")
    if os.path.exists(_meta_path):
        with open(_meta_path) as _fh:
            _ML_METADATA = json.load(_fh)
except Exception as _e:
    _HAVE_ML = False
    _ML_BUNDLE = None
    _ML_METADATA = {}

# ─── Load PSSM ───
PSSM = A.load_pssm()

# ─── Physicochemical scales (for visualization) ───
EISENBERG_SCALE = {
    'A': 0.62,'R':-2.53,'N':-0.78,'D':-0.90,'C': 0.29,'Q':-0.85,'E':-0.74,
    'G': 0.48,'H':-0.40,'I': 1.38,'L': 1.06,'K':-1.50,'M': 0.64,'F': 1.19,
    'P': 0.12,'S':-0.18,'T':-0.05,'W': 0.81,'Y': 0.26,'V': 1.08,
}

AA_COLORS = {
    'A':'#F4A460','V':'#F4A460','L':'#F4A460','I':'#F4A460','M':'#F4A460',  # hydrophobic
    'F':'#2E8B57','W':'#2E8B57','Y':'#2E8B57',                              # aromatic
    'K':'#DC143C','R':'#DC143C','H':'#DC143C',                              # positive
    'D':'#4169E1','E':'#4169E1',                                            # negative
    'S':'#6495ED','T':'#6495ED','N':'#6495ED','Q':'#6495ED','C':'#6495ED',  # polar
    'G':'#9370DB','P':'#9370DB',                                            # special
}


# ═════════════════════════════════════════════════════════════════════════════
# Helper functions
# ═════════════════════════════════════════════════════════════════════════════
def parse_fasta_text(text):
    """Parse FASTA text into [(id, seq), ...]."""
    entries = []
    name, buf = None, []
    for line in text.strip().split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                entries.append((name, "".join(buf).upper()))
            name = line[1:].split()[0] if len(line) > 1 else "seq"
            buf = []
        else:
            buf.append(re.sub(r"[^A-Za-z]", "", line))
    if name is not None:
        entries.append((name, "".join(buf).upper()))
    elif buf:
        entries.append(("sequence_1", "".join(buf).upper()))
    return entries


def parse_csv_input(df):
    """Auto-detect columns and return [(id, seq, uniprot), ...]."""
    cols = [c.lower().strip() for c in df.columns]
    col_map = {}
    for i, c in enumerate(cols):
        if any(k in c for k in ["seq"]):
            col_map["seq"] = df.columns[i]
        elif any(k in c for k in ["uniprot", "accession", "entry"]):
            col_map["uni"] = df.columns[i]
        elif any(k in c for k in ["protein", "name", "gene", "label"]):
            col_map["name"] = df.columns[i]
    entries = []
    for _, row in df.iterrows():
        name = str(row.get(col_map.get("name", ""), "")).strip()
        uni  = str(row.get(col_map.get("uni", ""), "")).strip()
        seq  = str(row.get(col_map.get("seq", ""), "")).strip()
        seq  = re.sub(r"[^A-Za-z]", "", seq).upper()
        if not name:
            name = uni if uni and uni != "nan" else f"protein_{len(entries)+1}"
        if uni == "nan":
            uni = ""
        entries.append((name, seq, uni))
    return entries


def fetch_uniprot(accession):
    """Fetch sequence from UniProt by accession ID. Returns (seq, error)."""
    try:
        import requests
        resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", timeout=10)
        if resp.status_code == 200:
            lines = resp.text.strip().split("\n")
            seq = "".join(l.strip() for l in lines if not l.startswith(">"))
            return seq.upper(), None
        return None, f"HTTP {resp.status_code}"
    except Exception as e:
        return None, str(e)


def search_uniprot(query, organism="Human", reviewed=True, max_results=5):
    """
    Search UniProt by protein name / gene name.
    Returns list of dicts: [{accession, name, gene, organism, length, sequence}, ...]
    """
    try:
        import requests
    except ImportError:
        return [], "requests library not installed"

    # Build the query
    q_parts = [f'(protein_name:"{query}" OR gene:"{query}")']
    if organism:
        q_parts.append(f"(organism_name:{organism})")
    if reviewed:
        q_parts.append("(reviewed:true)")
    q = " AND ".join(q_parts)

    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": q,
        "format": "json",
        "fields": "accession,protein_name,gene_names,organism_name,length,sequence",
        "size": max_results,
    }
    try:
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code != 200:
            return [], f"HTTP {resp.status_code}"
        data = resp.json()
        results = []
        for entry in data.get("results", []):
            acc = entry.get("primaryAccession", "")
            pname = entry.get("proteinDescription", {}).get(
                "recommendedName", {}).get("fullName", {}).get("value", "")
            if not pname:
                sub = entry.get("proteinDescription", {}).get("submissionNames", [])
                pname = sub[0].get("fullName", {}).get("value", "") if sub else ""
            genes = entry.get("genes", [])
            gene = genes[0].get("geneName", {}).get("value", "") if genes else ""
            org = entry.get("organism", {}).get("scientificName", "")
            length = entry.get("sequence", {}).get("length", 0)
            seq = entry.get("sequence", {}).get("value", "")
            results.append(dict(
                accession=acc, name=pname, gene=gene,
                organism=org, length=length, sequence=seq,
            ))
        return results, None
    except Exception as e:
        return [], str(e)


def resolve_protein_input(name_or_id):
    """
    Given a string that could be a UniProt ID, gene name, or protein name,
    try to resolve it to a sequence. Returns (display_name, sequence, details, error).
    """
    text = name_or_id.strip()
    if not text:
        return None, None, None, "empty input"

    # 1. Check if it looks like a UniProt accession (e.g. P15311, Q9Y2D5)
    if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$', text) or \
       re.match(r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$', text) or \
       re.match(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]-\d+$', text):
        seq, err = fetch_uniprot(text)
        if seq:
            return text, seq, {"accession": text}, None
        # Fall through to name search

    # 2. Search by name/gene
    results, err = search_uniprot(text, organism="Human", reviewed=True, max_results=5)
    if err:
        return text, None, None, err
    if not results:
        # Try without reviewed filter
        results, err = search_uniprot(text, organism="", reviewed=False, max_results=5)
    if not results:
        return text, None, None, f"no UniProt results for '{text}'"

    # Return the top hit (caller may present choices to the user)
    top = results[0]
    display = top["gene"] or top["name"] or top["accession"]
    return display, top["sequence"], results, None


def run_screen(proteins, rii_thr, ri_thr, use_ml, ml_thr, strict):
    """Run AKAP screen using the shared backend in akap_screen.py.

    This keeps Streamlit consistent with the CLI: ML v2, length-aware correction,
    background_risk, and proteomic confidence all come from A.screen_protein().
    """
    from types import SimpleNamespace

    # Make the app slider affect the backend high-confidence threshold.
    # Very-high remains stricter than high.
    if use_ml:
        A.ML_HIGH_THR = float(ml_thr)
        A.ML_VHIGH_THR = max(float(ml_thr) + 0.10, 0.90)

    args = SimpleNamespace(
        ri_only=False,
        rii_only=False,
        rii_thr=rii_thr,
        ri_thr=ri_thr,
        strict=strict,
    )

    ml_bundle = _ML_BUNDLE if (use_ml and _HAVE_ML) else None
    rows = []
    for pid, seq in proteins:
        rows.extend(A.screen_protein(pid, seq, args, PSSM, ml_bundle=ml_bundle))

    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df

    # Backward-compatible aliases for existing app plotting code.
    if "start" not in df.columns and "win_start" in df.columns:
        df["start"] = df["win_start"]
    if "end" not in df.columns and "win_end" in df.columns:
        df["end"] = df["win_end"]
    if "helix_turns" not in df.columns and "contiguous_hydro_turns" in df.columns:
        df["helix_turns"] = df["contiguous_hydro_turns"]

    return df.sort_values(["passes_proteomic_filter", "pssm_score"], ascending=[False, False]).reset_index(drop=True)


# ═════════════════════════════════════════════════════════════════════════════
# Visualization functions
# ═════════════════════════════════════════════════════════════════════════════
def plot_helical_wheel(seq, title="Helical Wheel"):
    """Plot a helical wheel diagram for a peptide."""
    n = len(seq)
    angle_step = 100  # degrees per residue for α-helix
    fig = go.Figure()
    # Draw residues
    for i, aa in enumerate(seq):
        angle_rad = math.radians(i * angle_step - 90)
        r = 1.0 + (i * 0.02)  # slight spiral
        x = r * math.cos(angle_rad)
        y = r * math.sin(angle_rad)
        color = AA_COLORS.get(aa, "#999999")
        hydro = EISENBERG_SCALE.get(aa, 0)
        fig.add_trace(go.Scatter(
            x=[x], y=[y], mode="markers+text",
            marker=dict(size=max(30 - i, 18), color=color, opacity=0.8,
                        line=dict(width=2, color="white")),
            text=f"<b>{aa}</b><sub>{i+1}</sub>",
            textposition="middle center",
            textfont=dict(size=11, color="white"),
            hovertext=f"{aa}{i+1} | H={hydro:.2f}",
            hoverinfo="text",
            showlegend=False,
        ))
        # Connect consecutive residues
        if i > 0:
            prev_rad = math.radians((i-1) * angle_step - 90)
            r_prev = 1.0 + ((i-1) * 0.02)
            x_prev = r_prev * math.cos(prev_rad)
            y_prev = r_prev * math.sin(prev_rad)
            fig.add_trace(go.Scatter(
                x=[x_prev, x], y=[y_prev, y], mode="lines",
                line=dict(color="#cccccc", width=1), showlegend=False, hoverinfo="skip"
            ))
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        xaxis=dict(visible=False, range=[-1.8, 1.8]),
        yaxis=dict(visible=False, range=[-1.8, 1.8], scaleanchor="x"),
        width=380, height=380, margin=dict(l=10, r=10, t=40, b=10),
        plot_bgcolor="white",
    )
    return fig


def plot_hydrophobicity_profile(seq, title="Hydrophobicity Profile"):
    """Plot residue-by-residue hydrophobicity with helix periodicity overlay."""
    h = [EISENBERG_SCALE.get(aa, 0) for aa in seq]
    x = list(range(1, len(seq)+1))
    colors = [AA_COLORS.get(aa, "#999") for aa in seq]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=x, y=h, marker_color=colors, name="Hydrophobicity",
        hovertext=[f"{seq[i]}{i+1}: {h[i]:.2f}" for i in range(len(seq))],
        hoverinfo="text",
    ))
    # Add helix periodicity envelope (3.6 residue)
    if len(seq) >= 4:
        env_x = np.linspace(0, len(seq)-1, 200)
        env_y = 0.8 * np.cos(2 * np.pi * env_x / 3.6)
        fig.add_trace(go.Scatter(
            x=env_x+1, y=env_y, mode="lines", name="α-helix period (3.6)",
            line=dict(color="red", width=1.5, dash="dash"), opacity=0.5,
        ))
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        xaxis_title="Position", yaxis_title="Eisenberg Hydrophobicity",
        height=280, margin=dict(l=50, r=20, t=40, b=40),
        legend=dict(orientation="h", y=1.12),
    )
    return fig


def plot_score_gauge(score, max_score=25, label="PSSM Score"):
    """A gauge chart for the PSSM score."""
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=score,
        title={"text": label, "font": {"size": 14}},
        gauge=dict(
            axis=dict(range=[0, max_score]),
            bar=dict(color="#2E86C1"),
            steps=[
                dict(range=[0, 7], color="#FADBD8"),
                dict(range=[7, 12], color="#FDEBD0"),
                dict(range=[12, 18], color="#D5F5E3"),
                dict(range=[18, max_score], color="#82E0AA"),
            ],
            threshold=dict(line=dict(color="red", width=3), thickness=0.8, value=7),
        ),
    ))
    fig.update_layout(height=200, margin=dict(l=20, r=20, t=40, b=10))
    return fig


def plot_protein_map(seq, hits_df, protein_name):
    """Plot a linear protein map showing hit locations."""
    fig = go.Figure()
    L = len(seq)
    # Protein backbone
    fig.add_trace(go.Scatter(
        x=[1, L], y=[0, 0], mode="lines",
        line=dict(color="#BDC3C7", width=12), showlegend=False, hoverinfo="skip",
    ))
    # Hit regions
    colors = {"RII": "#3498DB", "RI": "#E74C3C"}
    for _, h in hits_df.iterrows():
        fig.add_trace(go.Scatter(
            x=[h["start"], h["end"]], y=[0, 0], mode="lines",
            line=dict(color=colors.get(h["isoform"], "#888"), width=16),
            name=f'{h["isoform"]} ({h["start"]}-{h["end"]})',
            hovertext=f'{h["isoform"]} | {h["start"]}-{h["end"]} | score={h["pssm_score"]}',
            hoverinfo="text",
        ))
    fig.update_layout(
        title=dict(text=f"{protein_name} ({L} aa)", font=dict(size=13)),
        xaxis_title="Residue", yaxis=dict(visible=False, range=[-0.5, 0.5]),
        height=120, margin=dict(l=50, r=20, t=35, b=35),
        legend=dict(orientation="h", y=1.3),
    )
    return fig


# ═════════════════════════════════════════════════════════════════════════════
# Streamlit App
# ═════════════════════════════════════════════════════════════════════════════
def main():
    st.set_page_config(
        page_title="AKAP Domain Screener",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    if "resolved_proteins" not in st.session_state:
        st.session_state["resolved_proteins"] = []
    if "last_uniprot_queries" not in st.session_state:
        st.session_state["last_uniprot_queries"] = ""

    # ── Header ──
    st.markdown("""
    # 🧬 AKAP Domain Screener
    **Screen proteins for A-Kinase Anchoring Protein (AKAP) amphipathic-helix motifs**
    **with DDIP vs AKAP discrimination**

    Based on:
    - Burgers et al. (2015) *Biochemistry* 54, 11–21 — PSSM from THAHIT tool
    ([DOI: 10.1021/bi500721a](https://doi.org/10.1021/bi500721a))
    - Falcone & Scott (2025) *Biochem J* 482, 485–498 — DDIP/AKAP classification
    ([DOI: 10.1042/BCJ20253085](https://doi.org/10.1042/BCJ20253085))
    """)

    with st.expander("📖 Quick start: how to use this tool", expanded=False):
        st.markdown("""
**Step 1**: Enter your protein(s) — paste a sequence, upload a FASTA/CSV, or search UniProt by name.

**Step 2**: Click **Run AKAP Screen**.

**Step 3**: Look at the **proteomic_confidence** column in the results:

| You see | It means | What to do |
|---|---|---|
| 🟢 **Very High** | Strong AKAP candidate | ✅ Prioritize for experimental follow-up |
| 🟡 **High** | Likely AKAP | ✅ Good candidate for validation |
| 🔵 **Medium** | Moderate evidence | ❓ Review the details before deciding |
| ⚪ **Sensitive Only** | Weak match (PSSM only) | ⚠️ Probably NOT an AKAP — could be a transmembrane helix |
| 🔴 **Unlikely** | Biological red flags | ❌ Not an AKAP |

**Tip**: For large protein lists (>100 proteins), switch to **ML-prioritized proteomic confidence** mode in the sidebar — it filters out false positives automatically.
        """)

    # ── Sidebar: Parameters ──
    with st.sidebar:
        st.header("⚙️ Parameters")

        st.subheader("🎯 Screening Mode")
        screening_mode = st.radio(
            "Select mode:",
            ["Sensitive discovery", "ML-prioritized proteomic confidence"],
            help="Sensitive = more candidates, higher FPs. "
                 "ML-prioritized = stricter, uses ML + length-aware background risk for large protein lists.",
        )
        if screening_mode == "Sensitive discovery":
            st.info("💡 **Sensitive mode**: Shows ALL hits above the PSSM threshold. "
                    "Many will be false positives (transmembrane helices, random hydrophobic stretches). "
                    "Use the `proteomic_confidence` column to decide which hits to trust.")
        else:
            st.success("🔬 **Proteomic mode**: Shows only high/very_high confidence hits. "
                       "These have strong PSSM scores, ML support, and pass biological checks. "
                       "Best for screening large protein lists where you want reliable candidates.")

        st.subheader("PSSM Thresholds")
        rii_thr = st.slider("PKA-RIIα threshold", 3.0, 15.0, 7.0, 0.5,
                            help="Default 7.0 recovers 99.2% of known motifs")
        ri_thr = st.slider("PKA-RIα threshold", 5.0, 25.0, 12.0, 0.5,
                           help="Default 12.0 recovers 100% of known motifs")
        st.caption("💡 Lower thresholds = more hits but more false positives. "
                   "The defaults (7.0 / 12.0) are calibrated for sensitive discovery. "
                   "The proteomic confidence tiers apply stricter cutoffs automatically.")

        st.subheader("Options")
        strict = st.checkbox("Require literal regex match", value=False,
                             help="Only report hits that also match the consensus regex (lower recall)")

        use_ml = st.checkbox("Enable ML scoring", value=_HAVE_ML,
                             disabled=not _HAVE_ML,
                             help="Use the selected v2 ML model as the second-stage proteomic prioritization layer")
        ml_thr = 0.0
        if use_ml and _HAVE_ML:
            ml_thr = st.slider("ML high-confidence threshold", 0.0, 1.0, 0.80, 0.05)
            st.caption("💡 This sets the ML probability cutoff for 'high' confidence. "
                       "Hits below this are downgraded. Default 0.80 balances sensitivity and specificity.")

        st.divider()
        st.subheader("📊 About the tool")
        st.markdown("""
        **PSSM score**: Log-odds profile from 849 RII + 28 RI motifs. Higher = more AKAP-like.

        **ML probability**: second-stage proteomic prioritization model using interpretable
        sequence-derived features (PSSM score, hydrophobic moment, amphipathicity,
        helix propensity, anchor hydrophobicity, charge/composition features).

        **Length-aware background risk** accounts for the fact that long proteins have
        more sliding-window chances to produce random high-scoring PSSM hits.

        Validation metrics are loaded from `akap_ml_model_v2_metadata.json` when available;
        do not treat internal metrics as external biological validation.

        ---

        **🆕 DDIP vs AKAP classification**
        *(Falcone & Scott, Biochem J 2025)*

        | Class | Meaning |
        |---|---|
        | **AKAP** | Strong PKA anchor (≥5 helix turns, no negdets) |
        | **DDIP** | D/D domain interactor — binds d/d groove but NOT PKA |
        | **ambiguous** | Could be either — needs experimental validation |
        | **unlikely** | Charged residue (D/E) disrupts hydrophobic face |

        **Negative determinants**: Asp/Glu at hydrophobic anchor positions
        abolish PKA binding (demonstrated in OPA1 fungal forms & smAKAP S66D).
        """)

    # ── Input ──
    st.header("📥 Input")
    input_method = st.radio(
        "Choose input method:",
        ["Paste sequence(s)", "Upload FASTA file", "Upload CSV file",
         "Search by protein name / UniProt ID"],
        horizontal=True,
    )

    proteins = []  # list of (name, sequence)

    if input_method == "Paste sequence(s)":
        text = st.text_area(
            "Paste FASTA or raw sequence(s):",
            height=180,
            placeholder=">Protein1\nMAADSGRLH...\n>Protein2\nMKWVTFIS...\n\n(or just paste a raw sequence)",
        )
        if text.strip():
            if text.strip().startswith(">"):
                proteins = parse_fasta_text(text)
            else:
                seq = re.sub(r"[^A-Za-z]", "", text).upper()
                proteins = [("input_sequence", seq)]

    elif input_method == "Upload FASTA file":
        uploaded = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "faa", "txt"])
        if uploaded:
            text = uploaded.read().decode("utf-8")
            proteins = parse_fasta_text(text)
            st.success(f"Loaded {len(proteins)} sequence(s)")

    elif input_method == "Upload CSV file":
        uploaded = st.file_uploader("Upload CSV file", type=["csv", "tsv", "txt"])
        if uploaded:
            sep = "\t" if uploaded.name.endswith(".tsv") else ","
            df_in = pd.read_csv(uploaded, sep=sep, dtype=str).fillna("")
            st.write("Preview:", df_in.head())
            entries = parse_csv_input(df_in)
            # Separate: has sequence vs needs fetching
            have_seq = [(e[0], e[1]) for e in entries if len(e[1]) >= 20]
            need_resolve = [e for e in entries if len(e[1]) < 20]
            proteins = list(have_seq)

            if need_resolve:
                st.info(f"🔍 Resolving {len(need_resolve)} protein(s) via UniProt...")
                prog = st.progress(0)
                for idx, (name, _, uni) in enumerate(need_resolve):
                    # Try UniProt ID first, then name search
                    query = uni if uni else name
                    display, seq, details, err = resolve_protein_input(query)
                    if seq:
                        proteins.append((display or name, seq))
                        st.caption(f"  ✅ {query} → {display} ({len(seq)} aa)")
                    else:
                        st.warning(f"  ❌ Could not resolve '{query}': {err}")
                    prog.progress((idx+1) / len(need_resolve))
            st.success(f"{len(proteins)} protein(s) ready for screening")

    elif input_method == "Search by protein name / UniProt ID":
        st.markdown("""
        Enter **protein names**, **gene names**, or **UniProt accession IDs** — one per line.
        The app will search UniProt and fetch the sequences automatically.

        *Examples: `Ezrin`, `AKAP10`, `WAVE1`, `P15311`, `O43572`*
        """)

        query_text = st.text_area(
            "Enter protein/gene names or UniProt IDs:",
            placeholder="Ezrin\nAKAP10\nWASP family member 1\nP15311\nSPHKAP",
            height=150,
        )

        # Organism filter
        col_org1, col_org2 = st.columns([1, 1])
        with col_org1:
            organism = st.selectbox("Organism", ["Human", "Mouse", "Rat", "Any"], index=0)
        with col_org2:
            reviewed_only = st.checkbox("Reviewed (Swiss-Prot) only", value=True)

        if query_text.strip():
            queries = [q.strip() for q in query_text.strip().split("\n") if q.strip()]

            if st.button(f"🔍 Search UniProt ({len(queries)} queries)", use_container_width=True):
                st.divider()
                prog = st.progress(0)

                # Store search results in session state for user selection
                all_search_results = []
                for idx, query in enumerate(queries):
                    st.markdown(f"**Searching: `{query}`**")

                    # Check if it looks like a UniProt accession
                    is_accession = bool(
                        re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$', query) or
                        re.match(r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$', query)
                    )

                    if is_accession:
                        seq, err = fetch_uniprot(query)
                        if seq:
                            proteins.append((query, seq))
                            st.success(f"  ✅ {query}: fetched ({len(seq)} aa)")
                        else:
                            st.warning(f"  ❌ {query}: {err}")
                    else:
                        org = "" if organism == "Any" else organism
                        results, err = search_uniprot(query, organism=org,
                                                       reviewed=reviewed_only, max_results=5)
                        if err:
                            st.warning(f"  ❌ Search error: {err}")
                        elif not results:
                            # Retry without reviewed filter
                            results, _ = search_uniprot(query, organism=org,
                                                         reviewed=False, max_results=5)
                        if results:
                            # Show results as a table for transparency
                            res_df = pd.DataFrame([
                                {"": "✅" if i == 0 else "",
                                 "Accession": r["accession"],
                                 "Protein": r["name"][:50],
                                 "Gene": r["gene"],
                                 "Organism": r["organism"],
                                 "Length": r["length"]}
                                for i, r in enumerate(results)
                            ])
                            st.dataframe(res_df, use_container_width=True, hide_index=True)

                            # Use top hit
                            top = results[0]
                            display = top["gene"] or top["accession"]
                            proteins.append((f'{display}_{top["accession"]}', top["sequence"]))
                            st.success(
                                f"  → Using top hit: **{top['name'][:60]}** "
                                f"({top['accession']}, {top['gene']}, {top['length']} aa)")
                        else:
                            st.warning(f"  ❌ No UniProt results for '{query}'")

                    prog.progress((idx + 1) / len(queries))
                    time.sleep(0.3)  # Be polite to API

                if proteins:
                    st.session_state["resolved_proteins"] = list(proteins)
                    st.session_state["last_uniprot_queries"] = query_text
                    st.divider()
                    st.success(f"✅ {len(proteins)} protein(s) ready for screening. "
                              f"Click **Run AKAP Screen** below.")

    # Keep UniProt search results across Streamlit reruns.
    if input_method == "Search by protein name / UniProt ID" and not proteins:
        proteins = st.session_state.get("resolved_proteins", [])
        if proteins:
            st.info(f"Using {len(proteins)} resolved UniProt protein(s) from the previous search.")

    # ── Run screening ──
    if proteins:
        if st.button("🔍 Run AKAP Screen", type="primary", use_container_width=True):
            with st.spinner("Screening..."):
                t0 = time.time()
                results = run_screen(proteins, rii_thr, ri_thr, use_ml, ml_thr, strict)
                elapsed = time.time() - t0

            st.header("📊 Results")

            # Apply screening mode filter
            if len(results) > 0 and screening_mode == "ML-prioritized proteomic confidence":
                results = results[results.get("passes_proteomic_filter", pd.Series([False]*len(results))) == True].reset_index(drop=True)

            hit_proteins = set(results["protein"]) if len(results) > 0 else set()
            all_proteins = set(p[0] for p in proteins)
            no_hit = all_proteins - hit_proteins

            # Summary metrics
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Proteins screened", len(proteins))
            col2.metric("Proteins with hit", len(hit_proteins))
            col3.metric("Total motif hits", len(results))
            col4.metric("Time", f"{elapsed:.1f}s")

            # Confidence tier summary
            if len(results) > 0 and "proteomic_confidence" in results.columns:
                conf_counts = results["proteomic_confidence"].value_counts()
                tier_cols = st.columns(5)
                for i, tier in enumerate(["very_high", "high", "sensitive_only", "unlikely"]):
                    ct = conf_counts.get(tier, 0)
                    tier_cols[i].metric(tier.replace("_", " ").title(), ct)
                tier_cols[4].metric("Mode", screening_mode.split()[0])

                # ── User guide: what do these tiers mean? ──
                with st.expander("📖 What do these confidence tiers mean?", expanded=False):
                    st.markdown("""
| Tier | Meaning | Action |
|---|---|---|
| 🟢 **Very High** | Strong PSSM + ML support, low background noise | ✅ High-confidence AKAP candidate — prioritize for follow-up |
| 🟡 **High** | Good combined evidence (PSSM ≥ 12, ML ≥ 0.80) | ✅ Likely AKAP — consider for experimental validation |
| 🔵 **Medium** | Moderate ML support but below strict thresholds | ❓ Possible candidate — review manually before concluding |
| ⚪ **Sensitive Only** | Passed PSSM but failed ML / high background risk | ⚠️ Probably NOT an AKAP — could be a transmembrane helix or random match |
| 🔴 **Unlikely** | Charged residue (D/E) disrupts the hydrophobic face | ❌ Almost certainly NOT an AKAP |

**Key points for users:**
- **Very High + High** = results you can trust for proteomic-scale screening
- **Sensitive Only** = most of these are false positives; useful only if you have independent experimental evidence for that specific protein
- **ML unavailable** in filter reason = the ML model was not loaded; confidence is based on sequence rules only (less reliable)
- **Background risk = high** = this protein is long, so random hydrophobic stretches can match the PSSM pattern by chance
                    """)

                # Quick summary sentence
                n_high = conf_counts.get("very_high", 0) + conf_counts.get("high", 0)
                n_sens = conf_counts.get("sensitive_only", 0) + conf_counts.get("medium", 0)
                if n_high > 0:
                    st.success(f"🎯 **{n_high} high-confidence AKAP candidate(s)** found. "
                               f"{n_sens} additional sensitive-only hit(s).")
                elif n_sens > 0:
                    st.warning(f"⚠️ No high-confidence hits. {n_sens} sensitive-only hit(s) found — "
                               f"these are likely false positives unless you have prior evidence.")

            if len(results) > 0:
                # ── Results table ──
                st.subheader("🎯 Hit Table")

                # Format for display — proteomic confidence first
                display_cols = ["protein", "proteomic_confidence", "isoform", "classification",
                                "dual", "start", "end", "core", "pssm_score",
                                "ml_prob", "ml_confidence_tier", "passes_ml_filter",
                                "length_adjusted_score", "background_risk",
                                "amphipathic", "pI", "n_negdet", "helix_turns", "dd_class",
                                "filter_reason"]

                st.dataframe(
                    results[[c for c in display_cols if c in results.columns]],
                    use_container_width=True,
                    column_config={
                        "pssm_score": st.column_config.NumberColumn("PSSM Score", format="%.2f"),
                        "ml_prob": st.column_config.ProgressColumn("ML Prob", min_value=0, max_value=1, format="%.3f"),
                        "pI": st.column_config.NumberColumn("pI", format="%.2f"),
                        "helix_turns": st.column_config.NumberColumn("Helix Turns", format="%.1f"),
                        "length_adjusted_score": st.column_config.NumberColumn("Len-adjusted", format="%.2f"),
                        "passes_ml_filter": st.column_config.CheckboxColumn("ML Pass"),
                        "background_risk": st.column_config.TextColumn("Background Risk"),
                        "ml_confidence_tier": st.column_config.TextColumn("ML Tier"),
                        "dual": st.column_config.CheckboxColumn("Dual"),
                        "amphipathic": st.column_config.CheckboxColumn("Amphipathic"),
                        "proteomic_confidence": st.column_config.TextColumn("Confidence",
                            help="very_high/high = proteomic-grade, sensitive_only = candidate, unlikely = FP"),
                        "classification": st.column_config.TextColumn("Class",
                            help="AKAP=strong PKA anchor, DDIP=D/D domain interactor (no PKA), "
                                 "unlikely=charged residue disrupts hydrophobic face"),
                        "n_negdet": st.column_config.NumberColumn("NegDet",
                            help="Charged residues (D/E) at hydrophobic anchor positions — blocks PKA binding"),
                        "negdet_severity": st.column_config.TextColumn("Severity"),
                        "dd_class": st.column_config.TextColumn("d/d Class",
                            help="Predicted d/d domain partner: RIID2 (PKA-RII), RID2 (PKA-RI), DPY-30"),
                    },
                )

                # Download button
                csv_buf = results.to_csv(index=False)
                st.download_button(
                    "⬇️ Download results CSV", csv_buf, "akap_hits.csv", "text/csv",
                    use_container_width=True,
                )

                # ── Column guide ──
                with st.expander("📖 What do the columns mean?", expanded=False):
                    st.markdown("""
| Column | Meaning |
|---|---|
| **Confidence** | Overall verdict: very_high / high / medium / sensitive_only / unlikely |
| **PSSM Score** | How well the sequence matches known AKAP motifs (higher = better; ≥12 = strong) |
| **ML Prob** | Machine learning probability of being a real AKAP (0–1; ≥0.80 = good) |
| **ML Pass** | ✅ = ML supports this as an AKAP; ❌ = ML does not support it |
| **ML Tier** | ML-based confidence category |
| **Len-adjusted** | PSSM score corrected for protein length (longer proteins get more random matches) |
| **Background Risk** | How likely a random protein of this length would produce this PSSM score by chance |
| **Class** | AKAP (likely PKA anchor) / DDIP (binds d/d groove but NOT PKA) / unlikely (disrupted) |
| **NegDet** | Count of charged residues (D/E) at hydrophobic anchor positions — these block PKA binding |
| **d/d Class** | Which PKA regulatory subunit type: RIID2 (PKA-RIIα) or RID2 (PKA-RIα) |
| **Amphipathic** | ✅ = the helix has a proper hydrophobic + hydrophilic face pattern |
| **Helix Turns** | Estimated amphipathic helix length (AKAPs need ≥5 turns; DDIPs have 3–4) |
| **pI** | Isoelectric point of the motif core (most AKAPs are acidic, pI 3–6) |
| **Filter Reason** | Detailed explanation of why this hit got its confidence tier |
                    """)

                # ── Per-protein detail view ──
                st.subheader("🔬 Detailed View")
                sel_protein = st.selectbox("Select protein:", sorted(hit_proteins))
                sel_hits = results[results["protein"] == sel_protein]
                sel_seq = dict(proteins)[sel_protein]

                # Protein map
                st.plotly_chart(
                    plot_protein_map(sel_seq, sel_hits, sel_protein),
                    use_container_width=True,
                )

                # Per-hit visualizations
                for idx, (_, hit) in enumerate(sel_hits.iterrows()):
                    with st.expander(
                        f"**{hit['isoform']}** @ {hit['start']}–{hit['end']}  |  "
                        f"PSSM={hit['pssm_score']:.1f}"
                        + (f"  |  ML={hit['ml_prob']:.3f}" if pd.notna(hit.get('ml_prob', None)) else ""),
                        expanded=(idx == 0),
                    ):
                        core_seq = hit["core"]
                        win_seq = hit["window"]

                        c1, c2 = st.columns([1, 1])
                        with c1:
                            st.plotly_chart(
                                plot_helical_wheel(core_seq,
                                    f"Helical Wheel — {hit['isoform']} core"),
                                use_container_width=True,
                            )
                            st.caption("💡 In a true AKAP, hydrophobic residues (orange) cluster on one side "
                                       "and polar/charged residues (red/blue) on the other — this is the amphipathic pattern.")
                        with c2:
                            st.plotly_chart(
                                plot_hydrophobicity_profile(win_seq,
                                    f"Hydrophobicity — full window"),
                                use_container_width=True,
                            )
                            st.caption("💡 The red dashed line shows the α-helix periodicity (3.6 residues/turn). "
                                       "A real AKAP helix shows hydrophobic peaks aligning with this period.")

                        # Score gauge
                        max_s = 25 if hit["isoform"] == "RII" else 40
                        st.plotly_chart(
                            plot_score_gauge(hit["pssm_score"], max_s),
                            use_container_width=True,
                        )

                        # Details
                        detail_cols = st.columns(6)
                        detail_cols[0].metric("Classification", hit.get("classification", "—"))
                        detail_cols[1].metric("pI", f"{hit['pI']:.2f}")
                        detail_cols[2].metric("Helix turns", f"{hit.get('helix_turns', 0):.1f}")
                        detail_cols[3].metric("Neg. determinants", hit.get("n_negdet", 0))
                        detail_cols[4].metric("Canonical", "✅" if hit["canonical"] else "❌")
                        detail_cols[5].metric("d/d class", hit.get("dd_class", "—"))

                        # Interpretation hint
                        conf = hit.get("proteomic_confidence", "")
                        if conf == "very_high":
                            st.success("🟢 **Very high confidence** — this is a strong AKAP candidate with robust PSSM + ML support.")
                        elif conf == "high":
                            st.success("🟡 **High confidence** — good evidence for AKAP function. Consider experimental follow-up.")
                        elif conf == "medium":
                            st.info("🔵 **Medium confidence** — moderate support. Review the helical wheel and hydrophobicity profile before concluding.")
                        elif conf == "sensitive_only":
                            st.warning("⚪ **Sensitive only** — this hit passed the PSSM threshold but lacks ML support or has high background risk. "
                                       "Most sensitive-only hits in random proteins are false positives (e.g., transmembrane helices).")
                        elif conf == "unlikely":
                            st.error("🔴 **Unlikely** — a charged residue (Asp/Glu) sits on the hydrophobic face, which blocks PKA binding.")

                # ── Score distribution ──
                if len(results) > 3:
                    st.subheader("📈 Score Distributions")
                    fig_dist = px.histogram(
                        results, x="pssm_score", color="isoform", nbins=30,
                        color_discrete_map={"RII": "#3498DB", "RI": "#E74C3C"},
                        barmode="overlay", opacity=0.7,
                        labels={"pssm_score": "PSSM Score"},
                    )
                    fig_dist.update_layout(height=300)
                    st.plotly_chart(fig_dist, use_container_width=True)

            # ── Proteins without hits ──
            if no_hit:
                st.subheader("❌ Proteins without AKAP motifs")
                st.caption("These proteins have no sequence windows scoring above the PSSM threshold. "
                           "This means no amphipathic helix matching the AKAP pattern was found. "
                           "Note: some AKAPs bind PKA via non-canonical mechanisms that this tool cannot detect.")
                for p in sorted(no_hit):
                    st.write(f"  • {p}")

    # ── Footer ──
    st.divider()
    st.markdown("""
    <div style="text-align:center; color:#888; font-size:0.85em;">
    🧬 AKAP Domain Screener — kowith@ccs.tsukuba.ac.jp<br>
    PSSM + ML + DDIP/AKAP classification pipeline •
    <a href="https://doi.org/10.1021/bi500721a">Burgers et al. 2015</a> •
    <a href="https://doi.org/10.1042/BCJ20253085">Falcone & Scott 2025</a>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
