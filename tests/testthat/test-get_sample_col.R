source(file.path("..", "..", "R", "utils.R"))

test_that("get_sample_col prefers sampleID over patientID regardless of column order", {
    df <- data.frame(
        patientID = c("P1", "P1", "P2", "P2"),
        sampleID  = c("S1", "S2", "S3", "S4"),
        Group     = c("Control", "Treatment", "Control", "Treatment"),
        stringsAsFactors = FALSE
    )
    expect_equal(get_sample_col(df), "sampleID")
})

test_that("get_sample_col matches accepted names case-insensitively", {
    df <- data.frame(SAMPLE = c("a", "b"), Group = c("A", "B"))
    expect_equal(get_sample_col(df), "SAMPLE")
})

test_that("get_sample_col errors clearly when no accepted column exists, instead of guessing", {
    df <- data.frame(patientID = c("P1", "P2"), Group = c("A", "B"))
    expect_error(get_sample_col(df), "No sample ID column found")
})

test_that("check_sample_overlap reports mismatched samples instead of dropping them silently", {
    result <- check_sample_overlap(c("S1", "S2", "S3"), c("S1", "S2", "S4"))
    expect_equal(sort(result$common), c("S1", "S2"))
    expect_match(result$message, "S3")
    expect_match(result$message, "S4")
})

test_that("check_sample_overlap errors when metadata and counts share no samples", {
    expect_error(
        check_sample_overlap(c("S1", "S2"), c("X1", "X2")),
        "No matching sample IDs"
    )
})
