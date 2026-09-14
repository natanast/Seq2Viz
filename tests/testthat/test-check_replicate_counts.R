source(file.path("..", "..", "R", "utils.R"))

test_that("fewer than 2 samples in either level blocks with the observed counts", {
    meta <- data.frame(Group = c("control", "treatment", "treatment"))
    result <- check_replicate_counts(meta, "Group", "control", "treatment")

    expect_equal(result$status, "blocking")
    expect_match(result$message, "control = 1")
    expect_match(result$message, "treatment = 2")
})

test_that("exactly 2 in the smaller level runs but is flagged severe", {
    meta <- data.frame(Group = c("control", "control", "treatment", "treatment", "treatment"))
    result <- check_replicate_counts(meta, "Group", "control", "treatment")

    expect_equal(result$status, "severe")
    expect_match(result$message, "underpowered")
    expect_match(result$message, "control = 2")
    expect_match(result$message, "treatment = 3")
})

test_that("3 to 5 in the smaller level runs with a mild note", {
    meta <- data.frame(Group = c(rep("control", 3), rep("treatment", 4)))
    result <- check_replicate_counts(meta, "Group", "control", "treatment")

    expect_equal(result$status, "mild")
    expect_match(result$message, "large")
    expect_match(result$message, "control = 3")
    expect_match(result$message, "treatment = 4")
})

test_that("more than 5 in the smaller level is unflagged", {
    meta <- data.frame(Group = c(rep("control", 6), rep("treatment", 8)))
    result <- check_replicate_counts(meta, "Group", "control", "treatment")

    expect_equal(result$status, "ok")
    expect_match(result$message, "control = 6")
    expect_match(result$message, "treatment = 8")
})

test_that("severity is judged by the smaller of the two levels, not the total", {
    # 2 in control, 20 in treatment: total is plenty, but the comparison is still
    # bottlenecked by the smaller group.
    meta <- data.frame(Group = c(rep("control", 2), rep("treatment", 20)))
    result <- check_replicate_counts(meta, "Group", "control", "treatment")

    expect_equal(result$status, "severe")
})
