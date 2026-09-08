source(file.path("..", "..", "R", "utils.R"))

test_that("non-numeric columns are always categorical", {
    expect_true(is_categorical_covariate(c("A", "B", "A", "B")))
    expect_true(is_categorical_covariate(factor(c("A", "B"))))
})

test_that("a numerically-coded batch with few levels relative to sample count is categorical", {
    batch <- rep(c(1, 2, 3), 4) # 3 distinct values across 12 samples
    expect_true(is_categorical_covariate(batch))
})

test_that("a continuous covariate with many distinct values is kept numeric", {
    age <- c(34, 41, 29, 55, 62, 38, 47, 51, 33, 60) # 10 distinct values across 10 samples
    expect_false(is_categorical_covariate(age))
})

test_that("a binary 0/1 numeric covariate is categorical", {
    sex <- rep(c(0, 1), 10)
    expect_true(is_categorical_covariate(sex))
})

test_that("a small study where every sample has a distinct numeric value stays continuous", {
    rin <- c(7.2, 8.1, 6.9, 7.8) # 4 distinct values across 4 samples -- low absolute
    # count, but no repetition, so this must NOT be treated as categorical
    expect_false(is_categorical_covariate(rin))
})
