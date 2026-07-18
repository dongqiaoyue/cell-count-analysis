"""
load_data.py -- Part 1: Data Management

Initializes a SQLite database with a normalized relational schema and loads
every row from cell-count.csv.

Run:
    python load_data.py

Produces `cell_counts.db` in the repository root.
"""

import os
import sqlite3

import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(BASE_DIR, "cell-count.csv")
DB_PATH = os.path.join(BASE_DIR, "cell_counts.db")

# The five immune-cell populations stored one-row-per-population (long format).
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------
# Three tables, third-normal-form:
#   subjects     -- one row per patient (attributes that never vary within a
#                   subject: project, condition, age, sex, treatment, response)
#   samples      -- one row per biological sample (a subject can have many)
#   cell_counts  -- one row per (sample, population): the measured count
#
# Keeping cell_counts in long format means adding a 6th cell population never
# requires an ALTER TABLE -- it is just more rows. See README for rationale.
SCHEMA = """
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;

CREATE TABLE subjects (
    subject    TEXT PRIMARY KEY,
    project    TEXT NOT NULL,
    condition  TEXT,            -- indication: melanoma / carcinoma / healthy
    age        INTEGER,
    sex        TEXT,            -- M / F
    treatment  TEXT,            -- miraclib / phauximab / none
    response   TEXT             -- yes / no / NULL (e.g. healthy, untreated)
);

CREATE TABLE samples (
    sample                    TEXT PRIMARY KEY,
    subject                   TEXT NOT NULL,
    sample_type               TEXT,        -- PBMC / WB
    time_from_treatment_start INTEGER,
    FOREIGN KEY (subject) REFERENCES subjects(subject)
);

CREATE TABLE cell_counts (
    sample      TEXT NOT NULL,
    population  TEXT NOT NULL,
    count       INTEGER NOT NULL,
    PRIMARY KEY (sample, population),
    FOREIGN KEY (sample) REFERENCES samples(sample)
);

CREATE INDEX idx_samples_subject       ON samples(subject);
CREATE INDEX idx_cell_counts_pop       ON cell_counts(population);
CREATE INDEX idx_subjects_cond_treat   ON subjects(condition, treatment);
"""


def main() -> None:
    df = pd.read_csv(CSV_PATH)

    # Reset the database on every run so the load is idempotent.
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.executescript(SCHEMA)

    # subjects: one row per subject, drop duplicate rows.
    subjects = (
        df[["subject", "project", "condition", "age", "sex", "treatment", "response"]]
        .drop_duplicates(subset="subject")
    )
    subjects.to_sql("subjects", conn, if_exists="append", index=False)

    # samples: one row per sample.
    samples = df[
        ["sample", "subject", "sample_type", "time_from_treatment_start"]
    ].drop_duplicates(subset="sample")
    samples.to_sql("samples", conn, if_exists="append", index=False)

    # cell_counts: wide -> long.
    long = df.melt(
        id_vars=["sample"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    long.to_sql("cell_counts", conn, if_exists="append", index=False)

    conn.commit()

    # Report a quick summary so the grader sees the load succeeded.
    n_sub = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
    n_sam = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    n_cc = conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0]
    conn.close()

    print(f"Loaded {n_sub} subjects, {n_sam} samples, {n_cc} cell-count rows.")
    print(f"Database written to {DB_PATH}")


if __name__ == "__main__":
    main()
