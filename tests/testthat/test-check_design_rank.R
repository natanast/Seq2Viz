source(file.path("..", "..", "R", "utils.R"))

test_that("check_design_rank passes silently for a full-rank paired design", {
    meta <- data.frame(
        patientID = c("P1", "P1", "P2", "P2", "P3", "P3"),
        Age       = c(30, 30, 45, 45, 60, 60),
        Group     = c("control", "treatment", "control", "treatment", "control", "treatment"),
        stringsAsFactors = FALSE
    )
    meta$patientID <- as.factor(meta$patientID)
    meta$Group <- as.factor(meta$Group)
    # Age is genuinely continuous here across the study even though it repeats per patient;
    # patientID and Group are NOT confounded (every patient has both levels of Group).

    expect_silent(check_design_rank(meta, ~ patientID + Group, c("patientID", "Group")))
})

test_that("check_design_rank detects a covariate fully determined by another variable and names both", {
    # Age is constant within each patientID -- the exact scenario from the brief -- so
    # patientID and Age are collinear once patientID is in the design.
    meta <- data.frame(
        patientID = c("P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"),
        Age       = c(30, 30, 45, 45, 60, 60, 52, 52),
        Group     = c("control", "treatment", "control", "treatment", "control", "treatment", "control", "treatment"),
        stringsAsFactors = FALSE
    )
    meta$patientID <- as.factor(meta$patientID)
    meta$Group <- as.factor(meta$Group)

    err <- tryCatch(
        check_design_rank(meta, ~ patientID + Age + Group, c("patientID", "Age", "Group")),
        error = function(e) e
    )

    expect_s3_class(err, "error")
    expect_match(err$message, "patientID")
    expect_match(err$message, "Age")
    expect_match(err$message, "confounded")
})

test_that("check_design_rank reports over-parameterisation distinctly when there are too many terms for the sample size", {
    # 4 samples, but Batch alone has one level per sample -- more parameters than samples.
    meta <- data.frame(
        Batch = c("B1", "B2", "B3", "B4"),
        Group = c("control", "control", "treatment", "treatment"),
        stringsAsFactors = FALSE
    )
    meta$Batch <- as.factor(meta$Batch)
    meta$Group <- as.factor(meta$Group)

    err <- tryCatch(
        check_design_rank(meta, ~ Batch + Group, c("Batch", "Group")),
        error = function(e) e
    )

    expect_s3_class(err, "error")
    expect_match(err$message, "parameters")
    expect_match(err$message, "sample")
    expect_false(grepl("confounded", err$message))
})

test_that("is_constant_within_groups requires the grouping variable to actually group samples", {
    # Every value of `group` is unique, so it cannot meaningfully "determine" anything --
    # this must not be reported as a trivial confound.
    expect_false(is_constant_within_groups(c(1, 2, 3, 4), c("a", "b", "c", "d")))
})

test_that("find_confounded_pair ignores variables that vary independently of each other", {
    # Fully crossed 2 (patient) x 2 (batch) x 2 (group) design: no variable is constant
    # within any level of another, so nothing here should be reported as confounded.
    meta <- data.frame(
        patientID = rep(c("P1", "P2"), each = 4),
        Batch     = rep(c("B1", "B1", "B2", "B2"), times = 2),
        Group     = rep(c("control", "treatment"), times = 4),
        stringsAsFactors = FALSE
    )
    expect_null(find_confounded_pair(meta, c("patientID", "Batch", "Group")))
})
