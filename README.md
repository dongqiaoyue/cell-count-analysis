# Loblaw Bio — Immune Cell Population Analysis

A small analytics pipeline for Bob Loblaw's clinical-trial cell-count data. It
loads the data into a relational SQLite database, computes cell-population
frequencies, tests which populations predict response to *miraclib*, and serves
an interactive Streamlit dashboard.

## Quick start

```bash
make setup       # install dependencies (pip install -r requirements.txt)
make pipeline    # init DB + load data (Part 1) + generate all outputs (Parts 2-4)
make dashboard   # launch the interactive dashboard
```

`make pipeline` runs `python load_data.py` then `python analysis.py`, with no
manual steps. All generated tables and the boxplot land in `outputs/`.

## Repository layout

```
load_data.py     Part 1 — schema definition + CSV load -> cell_counts.db
analysis.py      Parts 2-4 — frequency table, statistics, subset analysis
dashboard.py     Streamlit dashboard (reads cell_counts.db)
cell-count.csv   Source data (input)
Makefile         setup / pipeline / dashboard targets
requirements.txt Python dependencies
outputs/         Generated tables and plots (committed for reference)
```

## Database schema

Third-normal-form design with three tables:

| Table           | Grain                            | Key columns                                                               |
| --------------- | -------------------------------- | ------------------------------------------------------------------------- |
| `subjects`    | one row per patient              | `subject` (PK), project, condition, age, sex, treatment, response       |
| `samples`     | one row per biological sample    | `sample` (PK), `subject` (FK), sample_type, time_from_treatment_start |
| `cell_counts` | one row per (sample, population) | (`sample`, `population`) PK, count                                    |

**Rationale.**

- **Separating subjects from samples** removes redundancy. A patient contributes
  many samples over time; storing age/sex/treatment/response once per subject
  (rather than repeating it on all 3+ sample rows) prevents update anomalies and
  keeps the data consistent. The source data confirms these attributes never
  vary within a subject.
- **`cell_counts` is stored long, not wide** — one row per cell population
  instead of five columns. This is the key design decision: adding a sixth
  population (or a new assay) becomes *inserting rows*, never an `ALTER TABLE`.
  It also makes the Part 2 frequency computation a simple `GROUP BY sample`, and
  lets analytics generalize over populations without hard-coding column names.
- **Indexes** on `samples.subject`, `cell_counts.population`, and
  `subjects(condition, treatment)` cover the filter/join patterns the analyses
  actually use (per-population lookups and cohort filtering).

**Scaling to hundreds of projects / thousands of samples.** The long, normalized
layout scales naturally:

- Add a `projects` table and demote `project` to a FK on `subjects` once project
  metadata (sponsor, protocol, dates) matters — the current single-column form is
  a deliberate simplification for one dataset.
- The composite PK / indexes keep cohort queries selective as row counts grow;
  the fact table (`cell_counts`) is the only one that grows with #populations ×
  #samples, and it stays narrow (3 columns).
- For heavy analytical workloads the same schema ports directly to a columnar
  warehouse (BigQuery/Snowflake/DuckDB): `cell_counts` becomes a partitioned
  fact table (partition by project/time), `subjects`/`samples` become dimensions.
  Long format is exactly what columnar engines and BI tools expect.

## What each part does

**Part 2 — frequency table** (`outputs/cell_frequencies.csv`): for every sample,
total cell count and each population's percentage. Columns: `sample`,
`total_count`, `population`, `count`, `percentage`. 52,500 rows (10,500 samples ×
5 populations).

**Part 3 — responders vs non-responders** (melanoma + miraclib + PBMC only):
per-population boxplots (`outputs/part3_boxplot.png`) and a two-sided
**Mann–Whitney U** test per population (`outputs/part3_statistics.csv`). A
non-parametric test is used because relative-frequency distributions are not
guaranteed normal and the test is robust to that.

Five populations are tested against the same cohort, so the raw p-values are
inflated by multiple comparisons. Reported alongside each is a
**Benjamini–Hochberg** adjusted q-value, which controls the false discovery
rate; only q-values are used to declare a hit.

Result (993 responder vs 975 non-responder samples):

| population           | median resp.    | median non-resp. | p-value           | q-value (BH)      | significant |
| -------------------- | --------------- | ---------------- | ----------------- | ----------------- | ----------- |
| b_cell               | 9.43            | 9.79             | 5.6e-02           | 1.4e-01           | no          |
| cd8_t_cell           | 24.73           | 24.60            | 6.4e-01           | 6.4e-01           | no          |
| **cd4_t_cell** | **30.22** | **29.66**  | **1.3e-02** | **6.7e-02** | no          |
| nk_cell              | 14.51           | 14.80            | 1.2e-01           | 2.0e-01           | no          |
| monocyte             | 19.61           | 19.94            | 1.6e-01           | 2.0e-01           | no          |

**CD4 T cells are the leading candidate biomarker**, but the evidence is
suggestive rather than conclusive. Responders carry a higher CD4 T-cell relative
frequency and the difference is nominally significant (p ≈ 0.013); after
correcting for the five tests it does not clear the 0.05 threshold (q ≈ 0.067).
B cells trend lower in responders (p ≈ 0.056, q ≈ 0.14) and are the second
candidate. The honest read: CD4 is worth prospective validation in an
independent cohort, and no population is confirmed on this data alone. The
effect sizes are also small relative to the spread — a difference in medians of
roughly half a percentage point — so a larger cohort, not just a different test,
is what would settle it.

**Part 4 — baseline subset** (melanoma / PBMC / miraclib / time=0):
`outputs/part4_baseline_subset.csv` plus breakdown counts in
`outputs/part4_counts.txt`. All three breakdowns are reported at the same grain
— one count per sample — so they are comparable and each sums to the subset
size. Here every subject contributes exactly one baseline sample (656 samples
from 656 subjects), so sample- and subject-level counts coincide.

- Samples per project — **prj1: 384, prj3: 272**
- Samples by response — **responders 331, non-responders 325**
- Samples by sex — **M 344, F 312**

**Starred question** — average B-cell count for melanoma males (all sample &
treatment types), responders, at time=0: **10206.15** (n=485 samples).

## Dashboard

`make dashboard` starts a local Streamlit server (default
`http://localhost:8501`) with three tabs mirroring Parts 2–4: a filterable
frequency table with CSV export, interactive responder-vs-non-responder boxplots
with the significance table, and the baseline subset breakdown. Because it reads
`cell_counts.db` live, it always reflects the loaded data.

> **Live dashboard:**
> https://cell-count-analysis-scljafbnjai7hhxvybsumd.streamlit.app/
>
> Hosted on Streamlit Community Cloud, deployed from this repository's `main`
> branch. `cell_counts.db` is committed, so the hosted app serves the same data
> without needing the pipeline to run there. It can also be run locally with
> `make dashboard`.
