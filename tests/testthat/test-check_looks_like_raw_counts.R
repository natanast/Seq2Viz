source(file.path("..", "..", "R", "utils.R"))

test_that("raw integer counts are not flagged", {
    mm <- matrix(c(0, 5, 10, 100, 0, 1, 250, 3), ncol = 2)
    result <- check_looks_like_raw_counts(mm)

    expect_true(result$looks_raw)
    expect_null(result$message)
})

test_that("TPM-like data (mostly fractional, many below 1) is flagged", {
    set.seed(1)
    mm <- matrix(c(0, 0.02, 0.15, 0.87, 3.4, 12.1, 0.05, 0.3, 45.6, 2.2), ncol = 2)
    result <- check_looks_like_raw_counts(mm)

    expect_false(result$looks_raw)
    expect_match(result$message, "not whole numbers")
    expect_match(result$message, "below 1")
})

test_that("a small amount of non-integer noise below the threshold is not flagged", {
    # 1 non-integer value out of 200 nonzero values (0.5%) -- below the 1% threshold.
    mm <- matrix(c(1:199, 7.3), ncol = 1)
    result <- check_looks_like_raw_counts(mm)

    expect_true(result$looks_raw)
})

test_that("an all-zero matrix is not flagged (nothing to judge)", {
    mm <- matrix(0, nrow = 5, ncol = 5)
    result <- check_looks_like_raw_counts(mm)

    expect_true(result$looks_raw)
})

test_that("normalised counts with many non-integers but few values below 1 are flagged for that reason only", {
    mm <- matrix(c(10.5, 20.25, 300.75, 4.1, 55.9, 0, 120.3, 8.8), ncol = 2)
    result <- check_looks_like_raw_counts(mm)

    expect_false(result$looks_raw)
    expect_match(result$message, "not whole numbers")
    expect_false(grepl("below 1", result$message))
})
