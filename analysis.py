"""
analysis.py -- Parts 2 and 3.

Reads cell_counts.db (created by load_data.py), computes every required
output table / plot, and writes them to ./outputs/.

Run:
    python analysis.py

Outputs:
    outputs/cell_frequencies.csv       Part 2 summary table
    outputs/part3_boxplot.png          Part 3 responder vs non-responder boxplots
    outputs/part3_statistics.csv       Part 3 Mann-Whitney U test per population
"""

import os
import sqlite3

import matplotlib

matplotlib.use("Agg")  # headless / Codespaces safe
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(BASE_DIR, "cell_counts.db")
OUT_DIR = os.path.join(BASE_DIR, "outputs")

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


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

        ax.boxplot([resp, nonr], labels=["responder", "non-responder"])
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

    fig.suptitle(
        "Melanoma / miraclib / PBMC: cell population frequency by response",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "part3_boxplot.png"), dpi=120)
    plt.close(fig)

    stats = pd.DataFrame(stats_rows)
    stats["p_value"] = stats["p_value"].map(lambda x: f"{x:.3e}")
    return stats


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    conn = connect()

    freq = build_frequency_table(conn)
    freq.to_csv(os.path.join(OUT_DIR, "cell_frequencies.csv"), index=False)

    stats = responder_comparison(conn, freq)
    stats.to_csv(os.path.join(OUT_DIR, "part3_statistics.csv"), index=False)

    conn.close()

    print("Part 2: cell_frequencies.csv written ({} rows).".format(len(freq)))
    print("Part 3: statistics ->")
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
