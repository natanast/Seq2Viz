# Seq2Viz — pre-submission work plan

## Context

Seq2Viz is an R Shiny app for RNA-seq differential expression analysis (DESeq2 backend)
with PCA, volcano and heatmap visualisation. We are preparing an Applications Note
submission to *Bioinformatics* (Oxford University Press).

The journal enforces two rules that shape this work:

1. Software must not require significant time investment to install, and must run
   under nearly all conditions on a wide range of machines.
2. All code and data must be publicly available in stable repositories.

The paper's novelty claim is being rewritten around a capability the app already has
but the manuscript did not describe: **GUI-based specification of arbitrary multi-factor
DESeq2 designs** — unlimited covariates, explicit reference/target contrast selection,
and a live preview of the resulting design formula. Work in this plan should protect and
strengthen that claim.

Please work through the phases in order. Ask before making architectural changes that
affect more than one module.

---

## Phase 1 — Correctness bugs (blocking)

### 1.1 `get_sample_col()` can silently select the wrong column

In the differential expression module:

```r
samp_col <- grep("sample|id", colnames(df), ignore.case = TRUE, value = TRUE)[1]
if (is.na(samp_col)) samp_col <- colnames(df)[1]
```

The `"id"` pattern also matches `patientID`. If a metadata file orders columns as
`patientID, sampleID, Group`, this returns `patientID`. The downstream
`match(sample_ids, as.character(meta[[samp_col]]))` then aligns counts to the wrong
samples, producing silently incorrect results with no error.

Required behaviour:

- Match against an explicit priority list of exact names (case-insensitive):
  `sampleID`, `sample_id`, `sample`, `Sample`.
- If none is found, raise a clear validation error naming the accepted column names.
- Remove the fallback to `colnames(df)[1]` entirely — guessing is what caused the bug.
- After matching, verify that the sample IDs in metadata and the count matrix column
  names intersect; if the overlap is empty or partial, report exactly which samples are
  missing from which file.

Add a regression test covering the `patientID`-before-`sampleID` column ordering.

### 1.2 Numeric covariates are coerced to factors

```r
for(col in design_cols) meta[[col]] <- as.factor(meta[[col]])
```

A continuous covariate (age, RIN, tumour purity) becomes one factor level per sample,
giving a rank-deficient model matrix.

Required behaviour: coerce to factor only when the column is non-numeric, or is numeric
with few distinct values relative to sample count. When a numeric column is kept
continuous, say so in the design formula preview so the user understands the model
being fitted.

### 1.3 Rank-deficient designs produce a cryptic error

A covariate confounded with the main factor triggers DESeq2's
"model matrix is not full rank", which is meaningless to a biologist.

Required behaviour: before constructing the `DESeqDataSet`, build the model matrix and
check its rank. If deficient, block the run and explain in plain language which variables
are confounded and that one of them should be removed from the design.

### 1.4 Minimum replicate check is too permissive

`need(nrow(meta) >= 2, ...)` allows a single sample per group, where DESeq2 cannot
estimate dispersion.

Required behaviour: require at least two samples in *each* level of the selected
contrast, and report the observed counts per level in the error message.

### 1.5 Silent `round()` on the count matrix

`mm <- round(mm)` will happily accept TPM, FPKM or already-normalised values and produce
meaningless output.

Required behaviour: detect input that does not look like raw integer counts (large
fraction of non-integers, or values below 1) and warn the user prominently before
proceeding. Do not block — some pipelines legitimately produce non-integer estimated
counts — but make the assumption visible.

---

## Phase 2 — Reproducibility export (highest-value new feature)

Add a **Download analysis report** control to the differential expression tab. The app
already holds everything needed; this is packaging, not new computation.

The export should be a plain-text or `.R` file containing:

- The full design formula as fitted, with covariates listed
- Main factor, reference level, target level
- Sample count per group after filtering, and the IDs of any samples dropped
- Shrinkage method (`ashr` via `lfcShrink`)
- Significance thresholds and visual parameters used for each generated plot
- Complete `sessionInfo()`

Ideally the output is a runnable script that reproduces the analysis outside the app from
the same input files. If full runnability is too large a change, a structured provenance
record is acceptable for now — flag which you implemented.

This underpins the reproducibility argument in the manuscript, so it should be
implemented cleanly rather than as an afterthought.

---

## Phase 3 — Installation and deployment

The current install path (manual `install.packages()` for ~16 CRAN packages, then
`BiocManager::install()` for 3 Bioconductor packages, then `setwd()` and `runApp()`)
conflicts directly with the journal's installation rule. All four items below are needed.

### 3.1 Dockerfile

Base it on a `rocker` image with R and Shiny preinstalled. Pin the R version. Install
system dependencies needed for the Bioconductor packages. Expose the Shiny port and set
the app as the entrypoint. Verify the image builds clean from scratch and the app is
reachable in a browser.

### 3.2 `renv.lock`

Run `renv::init()`, snapshot the full dependency tree, and commit the lockfile. Add a
short section to the README explaining `renv::restore()`.

### 3.3 Lower the R version requirement

The README requires R 4.5.3 — the latest patch release, which almost no shared cluster
or institutional installation will have. Unless a specific language feature requires it,
lower the floor to 4.3.0 or 4.4.0 and confirm the app runs there. Then make the stated
version consistent across README, Dockerfile and any in-app text.

### 3.4 One-step launch

`runApp()` from the repository root should work without `setwd()`. The README currently
instructs `runApp("R/main.R")` while an `app.R` exists at the root — resolve the
inconsistency in favour of the simplest path.

---

## Phase 4 — Interface and data handling

- The button labelled **"Download Results (Excel)"** writes CSV via `write.csv`. Either
  relabel it to CSV or implement genuine `.xlsx` export. Relabelling is preferred —
  it removes a dependency.
- Download buttons render before any analysis has run, offering empty files. Gate them
  behind a completed run.
- Provide plain-text (TSV/CSV) versions of everything in `example_data/`. Spreadsheet
  formats are a known hazard for gene identifiers, which Excel silently converts to
  dates. Keep the `.xlsx` files if you like, but the documented default should be TSV.

---

## Phase 5 — Repository hygiene for review

- Resolve or respond to the two open GitHub issues.
- Fix the README typo "Differential exprassion tab".
- Document the design-specification workflow in the README with a worked paired-design
  example, since this is now the paper's central claim.
- Document explicitly that DESeq2 is fitted on all samples while plots display only the
  two contrasted levels. This is a defensible choice — shared dispersion estimation
  across all samples is often preferable — but it is currently undocumented and a
  reviewer will notice.
- Add a minimal CI workflow that installs dependencies and confirms the app loads.
- Create a tagged release (`v1.0.0`) and connect the repository to Zenodo to mint a DOI.
  The manuscript has three `[Insert DOI]` placeholders that cannot be filled until this
  is done.

---

## Phase 6 — Optional, if time allows

- Deploy a public demo instance (shinyapps.io free tier is sufficient for review) with
  the example dataset preloaded. This answers most installation objections before they
  are raised.
- Verify rendering in Chrome, Firefox and Safari — the journal requires that web
  interfaces not be browser-specific.

---

## Working notes

- Do not weaken the design-specification flexibility while fixing the covariate handling;
  the multi-covariate `selectizeInput` and the live formula preview are the paper's
  central claim.
- Where a fix changes user-visible behaviour, note it so the manuscript's Methods section
  can be updated to match. The manuscript must describe what the code actually does.
