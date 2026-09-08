

# Accepted sample-ID column names, in priority order (matched case-insensitively).
# Deliberately does NOT match on a loose "id" substring: metadata files that also carry
# a patientID/subjectID column (common in paired designs) would otherwise have that
# column picked instead of the actual sample identifier.
SAMPLE_ID_COLUMN_CANDIDATES <- c("sampleID", "sample_id", "sample", "Sample")

# Identify the sample-ID column in a metadata data frame by exact name (case-insensitive),
# checked in priority order. Errors out rather than guessing when none of the accepted
# names is present.
get_sample_col <- function(df) {
    cols <- colnames(df)
    for (candidate in SAMPLE_ID_COLUMN_CANDIDATES) {
        hit <- cols[tolower(cols) == tolower(candidate)]
        if (length(hit) > 0) return(hit[1])
    }
    stop(
        "No sample ID column found in the metadata. Rename the column that identifies ",
        "samples to one of: ", paste(SAMPLE_ID_COLUMN_CANDIDATES, collapse = ", "),
        ". Columns present: ", paste(cols, collapse = ", "), "."
    )
}

# Verify that sample IDs found in metadata and in a second source (e.g. a count matrix's
# column names) actually overlap. Errors if there is no overlap at all (nothing to
# analyse); otherwise returns the common IDs plus a human-readable message describing any
# samples present in only one of the two sources, so callers can surface it to the user
# instead of silently dropping them.
check_sample_overlap <- function(meta_ids, other_ids, meta_label = "metadata", other_label = "count matrix") {
    meta_ids <- as.character(meta_ids)
    other_ids <- as.character(other_ids)

    common <- intersect(meta_ids, other_ids)
    missing_from_other <- setdiff(meta_ids, other_ids)
    missing_from_meta  <- setdiff(other_ids, meta_ids)

    if (length(common) == 0) {
        stop(
            "No matching sample IDs between ", meta_label, " and ", other_label, ". ",
            meta_label, " has: ", paste(meta_ids, collapse = ", "), ". ",
            other_label, " has: ", paste(other_ids, collapse = ", "), "."
        )
    }

    message <- NULL
    parts <- character(0)
    if (length(missing_from_other) > 0) {
        parts <- c(parts, paste0(
            length(missing_from_other), " sample(s) in ", meta_label, " missing from ",
            other_label, ": ", paste(missing_from_other, collapse = ", ")
        ))
    }
    if (length(missing_from_meta) > 0) {
        parts <- c(parts, paste0(
            length(missing_from_meta), " sample(s) in ", other_label, " missing from ",
            meta_label, ": ", paste(missing_from_meta, collapse = ", ")
        ))
    }
    if (length(parts) > 0) message <- paste(parts, collapse = " ")

    list(common = common, message = message)
}

# A numeric covariate is treated as categorical (coerced to a factor) only when it takes
# on few distinct values relative to the number of samples -- e.g. a batch numerically
# coded 1/2/3. Otherwise it is left as a continuous numeric covariate, since coercing a
# quantity like age, RIN, or tumour purity to a factor gives it one level per sample and
# produces a rank-deficient model matrix. Non-numeric columns are always categorical.
#
# Both thresholds must hold: at most NUMERIC_COVARIATE_MAX_LEVELS distinct values, AND
# those values covering at most NUMERIC_COVARIATE_MAX_FRACTION of the samples. The
# fraction check avoids misclassifying a genuinely continuous covariate in a small study
# where every sample happens to have a distinct value (e.g. 4 samples, 4 distinct RIN
# scores) as categorical just because the absolute count is low.
NUMERIC_COVARIATE_MAX_LEVELS <- 5
NUMERIC_COVARIATE_MAX_FRACTION <- 0.5

is_categorical_covariate <- function(x, max_levels = NUMERIC_COVARIATE_MAX_LEVELS,
                                      max_fraction = NUMERIC_COVARIATE_MAX_FRACTION) {
    if (!is.numeric(x)) return(TRUE)

    observed <- x[!is.na(x)]
    if (length(observed) == 0) return(TRUE)

    n_distinct <- length(unique(observed))
    n_distinct <= max_levels && (n_distinct / length(observed)) <= max_fraction
}

# TRUE if `x` is fully determined by `group` -- i.e. constant within every level of
# `group` -- which is how most real confounding between a covariate and a factor arises
# in practice (a numeric value that never varies within patient, or a batch that never
# varies within treatment group). Requires `group` to actually group more than one sample
# together; otherwise every `x` would be trivially "constant" within singleton groups,
# which would flag confounding between any variable and one that merely happens to be
# unique per sample (e.g. a nearly-continuous covariate) even though no real redundancy
# exists.
is_constant_within_groups <- function(x, group) {
    if (length(x) == 0) return(FALSE)
    if (length(unique(group)) >= length(group)) return(FALSE)

    groups <- split(x, group)
    all(vapply(groups, function(v) length(unique(v[!is.na(v)])) <= 1, logical(1)))
}

# Search all pairs of design columns for one being fully determined by the other. Returns
# c(determined_variable, determining_variable) for the first such pair found, or NULL if
# none of the pairs show this relationship.
find_confounded_pair <- function(meta, design_cols) {
    if (length(design_cols) < 2) return(NULL)

    combos <- utils::combn(design_cols, 2, simplify = FALSE)
    for (pair in combos) {
        a <- pair[1]; b <- pair[2]
        if (is_constant_within_groups(meta[[a]], meta[[b]])) return(c(a, b))
        if (is_constant_within_groups(meta[[b]], meta[[a]])) return(c(b, a))
    }
    NULL
}

# Checks a design matrix for rank deficiency before DESeq2 ever sees it, so a confounded
# or over-parameterised design produces a plain-language message instead of DESeq2's raw
# "model matrix is not full rank" error. Two distinct failure modes are reported
# separately, since they call for different fixes:
#   - Over-parameterisation: the formula has at least as many parameters as samples, e.g.
#     too many factor levels for the amount of data. Fix: simplify the design or collect
#     more samples.
#   - Confounding: two variables are collinear -- most commonly because one is fully
#     determined by the other, e.g. a numeric covariate that is constant within every
#     level of a factor. Fix: drop one of the confounded variables.
# Does nothing (returns invisibly) when the design is full rank.
check_design_rank <- function(meta, design_formula, design_cols) {
    mm <- stats::model.matrix(design_formula, data = as.data.frame(meta))
    n <- nrow(mm)
    p <- ncol(mm)

    if (p >= n) {
        stop(
            "The design '", paste(deparse(design_formula), collapse = " "), "' has ", p,
            " parameters but only ", n, " sample(s). There are too many factor levels ",
            "(or covariates) for the amount of data available to fit. Remove or simplify ",
            "one of: ", paste(design_cols, collapse = ", "), "."
        )
    }

    rank <- qr(mm)$rank
    if (rank == p) return(invisible(NULL))

    confounded <- find_confounded_pair(meta, design_cols)
    if (!is.null(confounded)) {
        stop(
            "'", confounded[1], "' and '", confounded[2], "' are confounded: ",
            confounded[1], " is fully determined by ", confounded[2],
            " (every sample with the same ", confounded[2], " has the same ", confounded[1],
            "), so the model cannot separate their effects. Remove one of them from the design."
        )
    }

    stop(
        "The design matrix is not full rank (rank ", rank, " of ", p, " parameters). ",
        "Some combination of variables in the design is redundant. Try removing one of: ",
        paste(design_cols, collapse = ", "), "."
    )
}
