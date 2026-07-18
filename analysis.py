"""
analysis.py -- Part 2.

Reads cell_counts.db (created by load_data.py), computes every required
output table / plot, and writes them to ./outputs/.

Run:
    python analysis.py

Outputs:
    outputs/cell_frequencies.csv       Part 2 summary table
"""

import os
import sqlite3

import pandas as pd

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


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    conn = connect()

    freq = build_frequency_table(conn)
    freq.to_csv(os.path.join(OUT_DIR, "cell_frequencies.csv"), index=False)

    conn.close()

    print("Part 2: cell_frequencies.csv written ({} rows).".format(len(freq)))


if __name__ == "__main__":
    main()
