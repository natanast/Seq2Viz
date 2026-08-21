

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
