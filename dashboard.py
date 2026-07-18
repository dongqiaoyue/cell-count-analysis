"""
dashboard.py -- interactive Streamlit dashboard for Bob's analysis.

Run:
    streamlit run dashboard.py     (or `make dashboard`)

Reads directly from cell_counts.db, so it always reflects the loaded data.
Run `python load_data.py` first if the database does not exist yet.
"""

import os
import sqlite3

import pandas as pd
import plotly.express as px
import streamlit as st
from scipy.stats import false_discovery_control, mannwhitneyu

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(BASE_DIR, "cell_counts.db")
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

st.set_page_config(page_title="Loblaw Bio -- Immune Cell Dashboard", layout="wide")


@st.cache_data
def load_tables():
    conn = sqlite3.connect(DB_PATH)
    cc = pd.read_sql("SELECT * FROM cell_counts", conn)
    samples = pd.read_sql("SELECT * FROM samples", conn)
    subjects = pd.read_sql("SELECT * FROM subjects", conn)
    conn.close()

    totals = cc.groupby("sample")["count"].sum().rename("total_count")
    freq = cc.join(totals, on="sample")
    freq["percentage"] = freq["count"] / freq["total_count"] * 100
    full = (
        freq.merge(samples, on="sample")
        .merge(subjects, on="subject")
    )
    return freq, full


if not os.path.exists(DB_PATH):
    st.error("cell_counts.db not found. Run `python load_data.py` first.")
    st.stop()

freq, full = load_tables()

st.title("Loblaw Bio — Immune Cell Population Dashboard")
st.caption("Interactive view of the cell-count clinical trial data.")

tab2, tab3, tab4 = st.tabs(
    ["Part 2 · Frequency table", "Part 3 · Responders", "Part 4 · Baseline subset"]
)

# --------------------------------------------------------------------------
# Part 2
# --------------------------------------------------------------------------
with tab2:
    st.subheader("Relative frequency of each cell population per sample")
    samples_sel = st.multiselect(
        "Filter by sample (leave empty for all)",
        sorted(freq["sample"].unique()),
    )
    table = freq[["sample", "total_count", "population", "count", "percentage"]].copy()
    if samples_sel:
        table = table[table["sample"].isin(samples_sel)]
    table["percentage"] = table["percentage"].round(2)
    st.dataframe(table, use_container_width=True, height=430)
    st.download_button(
        "Download CSV",
        table.to_csv(index=False),
        "cell_frequencies.csv",
        "text/csv",
    )

# --------------------------------------------------------------------------
# Part 3
# --------------------------------------------------------------------------
with tab3:
    st.subheader("Responders vs non-responders — melanoma / miraclib / PBMC")
    d = full[
        (full["condition"] == "melanoma")
        & (full["treatment"] == "miraclib")
        & (full["sample_type"] == "PBMC")
    ].copy()
    d["response"] = d["response"].map({"yes": "responder", "no": "non-responder"})

    fig = px.box(
        d,
        x="population",
        y="percentage",
        color="response",
        points="outliers",
        labels={"percentage": "relative frequency (%)"},
        category_orders={"population": POPULATIONS},
    )
    st.plotly_chart(fig, use_container_width=True)

    rows = []
    for pop in POPULATIONS:
        r = d[(d["population"] == pop) & (d["response"] == "responder")]["percentage"]
        n = d[(d["population"] == pop) & (d["response"] == "non-responder")]["percentage"]
        u, p = mannwhitneyu(r, n, alternative="two-sided")
        rows.append(
            {
                "population": pop,
                "median_responder": round(r.median(), 2),
                "median_nonresponder": round(n.median(), 2),
                "p_value": p,
            }
        )

    # Five populations tested on one cohort -> correct for multiple testing.
    q = false_discovery_control([r["p_value"] for r in rows], method="bh")
    for row, qv in zip(rows, q):
        row["q_value (BH)"] = f"{qv:.3e}"
        row["significant (FDR<0.05)"] = "✅" if qv < 0.05 else ""
        row["p_value"] = f"{row['p_value']:.3e}"

    st.markdown("**Mann–Whitney U test (two-sided), Benjamini–Hochberg corrected**")
    st.dataframe(pd.DataFrame(rows), use_container_width=True)
    st.caption(
        "Five populations are tested on the same cohort, so raw p-values are "
        "inflated. Only populations flagged ✅ survive false-discovery-rate "
        "correction and qualify as candidate predictors of miraclib response."
    )

# --------------------------------------------------------------------------
# Part 4
# --------------------------------------------------------------------------
with tab4:
    st.subheader("Melanoma PBMC baseline (t=0), miraclib-treated")
    sub = full[
        (full["condition"] == "melanoma")
        & (full["treatment"] == "miraclib")
        & (full["sample_type"] == "PBMC")
        & (full["time_from_treatment_start"] == 0)
    ]
    # All three breakdowns at one grain (per sample) so they stay comparable.
    samp = sub.drop_duplicates("sample")

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**Samples per project**")
        st.dataframe(samp["project"].value_counts().rename("samples"))
    with c2:
        st.markdown("**Samples by response**")
        st.dataframe(samp["response"].value_counts().rename("samples"))
    with c3:
        st.markdown("**Samples by sex**")
        st.dataframe(samp["sex"].value_counts().rename("samples"))
    st.caption(
        f"{len(samp)} baseline samples from {samp['subject'].nunique()} subjects "
        "— one baseline sample per subject, so sample and subject counts coincide."
    )

    starred = full[
        (full["condition"] == "melanoma")
        & (full["sex"] == "M")
        & (full["response"] == "yes")
        & (full["time_from_treatment_start"] == 0)
        & (full["population"] == "b_cell")
    ]
    st.metric(
        "Avg B-cell count — melanoma males, responders, t=0 (all sample/treatment types)",
        f"{starred['count'].mean():.2f}",
    )
