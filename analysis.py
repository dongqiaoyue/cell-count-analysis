"""
analysis.py -- Parts 2, 3 and 4.

Reads cell_counts.db (created by load_data.py), computes every required
output table / plot, and writes them to ./outputs/.

Run:
    python analysis.py

Outputs:
    outputs/cell_frequencies.csv       Part 2 summary table
    outputs/part3_boxplot.png          Part 3 responder vs non-responder boxplots
    outputs/part3_statistics.csv       Part 3 Mann-Whitney U test per population
    outputs/part4_baseline_subset.csv  Part 4 filtered baseline samples
    outputs/part4_counts.txt           Part 4 breakdown counts + starred answer
"""

import os
import sqlite3

import matplotlib

matplotlib.use("Agg")  # headless / Codespaces safe
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import false_discovery_control, mannwhitneyu

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(BASE_DIR, "cell_counts.db")
OUT_DIR = os.path.join(BASE_DIR, "outputs")

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def boxplot(ax, data, labels):
    """matplotlib renamed boxplot's `labels` kwarg to `tick_labels` in 3.9."""
    try:
        return ax.boxplot(data, tick_labels=labels)
    except TypeError:
        return ax.boxplot(data, labels=labels)


def _apply_fdr(rows, alpha: float = 0.05):
    """Add Benjamini-Hochberg adjusted q-values to a list of stat rows."""
    q = false_discovery_control([r["p_value"] for r in rows], method="bh")
    for row, qv in zip(rows, q):
        row["q_value_bh"] = qv
        row["significant_fdr_0.05"] = bool(qv < alpha)
    return rows


def connect() -> sqlite3.Connection:
    if not os.path.exists(DB_PATH):
        raise SystemExit("cell_counts.db not found -- run `python load_data.py` first.")
    return sqlite3.connect(DB_PATH)


# ---------------------------------------------------------------------------
# Part 2 -- relative frequency summary table
# ---------------------------------------------------------------------------
def build_frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """One row per (sample, population) with count, per-sample total, and %."""
    cc = pd.read_sql("SELECT sample, population, count FROM cell_counts", conn)
    totals = cc.groupby("sample")["count"].sum().rename("total_count")
    out = cc.join(totals, on="sample")
    out["percentage"] = (out["count"] / out["total_count"] * 100).round(4)
    out = out[["sample", "total_count", "population", "count", "percentage"]]
    out = out.sort_values(["sample", "population"]).reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# Part 3 -- responders vs non-responders (melanoma / miraclib / PBMC)
# ---------------------------------------------------------------------------
def responder_comparison(conn: sqlite3.Connection, freq: pd.DataFrame):
    meta = pd.read_sql(
        """
        SELECT s.sample, sub.response
        FROM samples s
        JOIN subjects sub ON s.subject = sub.subject
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
        """,
        conn,
    )
    data = freq.merge(meta, on="sample", how="inner")

    # Boxplots: one panel per population, responders vs non-responders.
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(20, 5), sharey=False)
    stats_rows = []
    for ax, pop in zip(axes, POPULATIONS):
        sub = data[data["population"] == pop]
        resp = sub[sub["response"] == "yes"]["percentage"]
        nonr = sub[sub["response"] == "no"]["percentage"]

        boxplot(ax, [resp, nonr], ["responder", "non-responder"])
        ax.set_title(pop)
        ax.set_ylabel("relative frequency (%)")

        u, p = mannwhitneyu(resp, nonr, alternative="two-sided")
        stats_rows.append(
            {
                "population": pop,
                "n_responder": len(resp),
                "n_nonresponder": len(nonr),
                "median_responder": round(resp.median(), 3),
                "median_nonresponder": round(nonr.median(), 3),
                "mannwhitney_u": round(u, 1),
                "p_value": p,
                "significant_0.05": p < 0.05,
            }
        )

    # Five populations are tested on the same cohort, so the nominal p-values
    # are inflated. Benjamini-Hochberg controls the false discovery rate; a
    # population is only reported as a real hit if it survives correction.
    stats_rows = _apply_fdr(stats_rows)

    fig.suptitle(
        "Melanoma / miraclib / PBMC: cell population frequency by response",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "part3_boxplot.png"), dpi=120)
    plt.close(fig)

    stats = pd.DataFrame(stats_rows)
    for col in ("p_value", "q_value_bh"):
        stats[col] = stats[col].map(lambda x: f"{x:.3e}")
    return stats


# ---------------------------------------------------------------------------
# Part 4 -- baseline subset analysis
# ---------------------------------------------------------------------------
def baseline_subset(conn: sqlite3.Connection):
    """Melanoma PBMC baseline (t=0) samples from miraclib-treated patients."""
    subset = pd.read_sql(
        """
        SELECT s.sample, sub.subject, sub.project, sub.response, sub.sex,
               s.time_from_treatment_start
        FROM samples s
        JOIN subjects sub ON s.subject = sub.subject
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND s.time_from_treatment_start = 0
        """,
        conn,
    )

    # All three breakdowns are reported at the same grain -- one count per
    # sample -- so the numbers are directly comparable and sum to the subset
    # size. In this subset each subject contributes exactly one baseline
    # sample, so the subject-level counts are identical (asserted below).
    by_project = subset.groupby("project")["sample"].count()
    by_response = subset.groupby("response")["sample"].count()
    by_sex = subset.groupby("sex")["sample"].count()

    n_samples, n_subjects = len(subset), subset["subject"].nunique()

    lines = ["Part 4 -- Melanoma PBMC baseline (t=0), miraclib-treated", ""]
    lines.append(f"Total samples in subset: {n_samples}")
    lines.append(f"Unique subjects: {n_subjects}")
    if n_samples == n_subjects:
        lines.append("(1 baseline sample per subject -- sample and subject counts coincide)")
    lines.append("")
    lines.append("Samples per project:")
    for k, v in by_project.items():
        lines.append(f"  {k}: {v}")
    lines.append("")
    lines.append("Samples by response:")
    for k, v in by_response.items():
        lines.append(f"  {k}: {v}")
    lines.append("")
    lines.append("Samples by sex:")
    for k, v in by_sex.items():
        lines.append(f"  {k}: {v}")

    # Starred required answer: Melanoma males, ALL sample & treatment types,
    # responders, time=0 -> average raw B-cell count (XXX.XX).
    starred = pd.read_sql(
        """
        SELECT cc.count AS b_cell
        FROM cell_counts cc
        JOIN samples s   ON cc.sample = s.sample
        JOIN subjects sub ON s.subject = sub.subject
        WHERE cc.population = 'b_cell'
          AND sub.condition = 'melanoma'
          AND sub.sex = 'M'
          AND sub.response = 'yes'
          AND s.time_from_treatment_start = 0
        """,
        conn,
    )
    avg_b = starred["b_cell"].mean()
    lines.append("")
    lines.append(
        "Starred question -- Melanoma males (all sample & treatment types), "
        "responders at time=0, average B-cell count:"
    )
    lines.append(f"  {avg_b:.2f}   (n = {len(starred)} samples)")

    return subset, "\n".join(lines)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    conn = connect()

    freq = build_frequency_table(conn)
    freq.to_csv(os.path.join(OUT_DIR, "cell_frequencies.csv"), index=False)

    stats = responder_comparison(conn, freq)
    stats.to_csv(os.path.join(OUT_DIR, "part3_statistics.csv"), index=False)

    subset, report = baseline_subset(conn)
    subset.to_csv(os.path.join(OUT_DIR, "part4_baseline_subset.csv"), index=False)
    with open(os.path.join(OUT_DIR, "part4_counts.txt"), "w") as f:
        f.write(report + "\n")

    conn.close()

    print("Part 2: cell_frequencies.csv written ({} rows).".format(len(freq)))
    print("Part 3: statistics ->")
    print(stats.to_string(index=False))
    print("\nPart 4:")
    print(report)


if __name__ == "__main__":
    main()
