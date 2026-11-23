test_that("detectionP, detectionP.neg, mask, eval_detP_cutoffs, correct_dye_bias, and beta_values work", {
  k450 = testthat::test_path("../../data/5640269011_R01C01")
  stopifnot(file.exists(paste0(k450, "_Grn.idat")))

  raw = read_idats(c(k450, k450), quiet = TRUE)

  # detectionP (linear scale)
  raw_p = detectionP(raw)
  expect_true("detP" %in% names(raw_p))
  expect_equal(dim(raw_p$detP), dim(raw$M))
  det_vals = raw_p$detP[!is.na(raw_p$detP)]
  expect_true(all(det_vals >= 0 & det_vals <= 1))

  # eval_detP_cutoffs with single male and single female (tests drop = FALSE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(eval_detP_cutoffs(raw_p, males = 1L, females = 2L))

  # detectionP.neg (log10 scale)
  raw_neg = detectionP.neg(raw)
  expect_equal(dim(raw_neg$detP), dim(raw$M))
  neg_vals = raw_neg$detP[!is.na(raw_neg$detP)]
  expect_true(all(neg_vals <= 0))

  # mask()
  raw_masked = mask(raw_p, 0.01)
  high_p = which(raw_p$detP[, 1] > 0.01)
  if (length(high_p) > 0) {
    expect_true(all(is.na(raw_masked$M[high_p, 1])))
    expect_true(all(is.na(raw_masked$U[high_p, 1])))
  }

  # correct_dye_bias() and beta_values()
  raw_corr = correct_dye_bias(raw)
  beta = beta_values(raw_corr)
  expect_true(is.matrix(beta))
  expect_equal(dim(beta), dim(raw$M))
  expect_equal(rownames(beta), as.character(raw$manifest$ilmn_id))
  expect_equal(colnames(beta), raw$meta$sample_id)
  b_vals = beta[!is.na(beta)]
  expect_true(all(b_vals > 0 & b_vals < 1))
})

test_that("check_sex and predict_sex infer sex and flag outliers", {
  k450 = testthat::test_path("../../data/5640269011_R01C01")
  stopifnot(file.exists(paste0(k450, "_Grn.idat")))

  raw = read_idats(k450, quiet = TRUE)
  xy = check_sex(raw)
  expect_named(xy, c("X", "Y"))
  expect_length(xy$X, 1)
  expect_length(xy$Y, 1)

  X = c(0.55, 0.57, 1.02, 1.05, 0.56, 1.03, 0.98)
  Y = c(1.10, 1.08, 0.15, 0.12, 1.12, 0.14, 0.95)
  pred = predict_sex(X, Y, male = 1:2, female = 3:4)

  expect_s3_class(pred, "factor")
  expect_equal(levels(pred), c("m", "f"))
  expect_equal(as.character(pred[1:6]), c("m", "m", "f", "f", "m", "f"))
  expect_true(is.na(pred[7]))
})

test_that("call_genotypes (learn = FALSE/TRUE), snp_outliers, and SNP donor agreement work", {
  set.seed(123)
  n_snps = 65

  draw_donor = function(n) {
    gt = sample(1:3, n, replace = TRUE)
    means = c(0.05, 0.50, 0.95)[gt]
    vals = stats::rnorm(n, mean = means, sd = 0.02)
    pmin(pmax(vals, 0.005), 0.995)
  }

  d1 = draw_donor(n_snps)
  d2 = draw_donor(n_snps)
  d3 = draw_donor(n_snps)
  contaminated = stats::runif(n_snps, min = 0.15, max = 0.85)

  snps = cbind(
    S1 = pmin(pmax(d1 + stats::rnorm(n_snps, 0, 0.005), 0.005), 0.995),
    S2 = pmin(pmax(d1 + stats::rnorm(n_snps, 0, 0.005), 0.005), 0.995),
    S3 = pmin(pmax(d2 + stats::rnorm(n_snps, 0, 0.005), 0.005), 0.995),
    S4 = pmin(pmax(d2 + stats::rnorm(n_snps, 0, 0.005), 0.005), 0.995),
    S5 = d3,
    S6 = contaminated
  )
  snps[1, 2] = NA_real_

  # Predefined parameters (learn = FALSE) and EM algorithm (learn = TRUE)
  gt_fixed = call_genotypes(snps, learn = FALSE)
  gt_em    = call_genotypes(snps, learn = TRUE, maxiter = 20)

  expect_null(gt_fixed$loglik)
  expect_true(length(gt_em$loglik) >= 2)
  expect_true(is.na(gt_em$outliers[1, 2]))
  expect_true(is.na(gt_em$gamma[[1]][1, 2]))

  gamma_sum = gt_fixed$gamma[[1]] + gt_fixed$gamma[[2]] + gt_fixed$gamma[[3]]
  expect_equal(gamma_sum[!is.na(gamma_sum)], rep(1, sum(!is.na(gamma_sum))), tolerance = 1e-8)

  # snp_outliers() flags contaminated sample S6
  out = snp_outliers(gt_fixed)
  expect_length(out, 6)
  expect_true(all(out[1:5] < -4))
  expect_true(out[6] > -4)

  # enumerate_sample_donors() and check_snp_agreement() on clean samples 1:5
  gt_clean = call_genotypes(snps[, 1:5], learn = FALSE)
  donors = enumerate_sample_donors(gt_clean)
  expect_equal(donors[1], donors[2])
  expect_equal(donors[3], donors[4])
  expect_false(donors[1] == donors[3])
  expect_false(donors[5] %in% donors[1:4])

  # No conflicts when donor labels match genotypes
  expect_null(check_snp_agreement(
    gt_clean,
    donor_ids = c("D1", "D1", "D2", "D2", "D3"),
    sample_ids = paste0("S", 1:5)
  ))

  # Conflicts detected with both character and integer sample_ids
  conflicts_chr = check_snp_agreement(
    gt_clean,
    donor_ids = c("D1", "D2", "D2", "D1", "D3"),
    sample_ids = paste0("S", 1:5)
  )
  conflicts_int = check_snp_agreement(
    gt_clean,
    donor_ids = c("D1", "D2", "D2", "D1", "D3"),
    sample_ids = 101:105
  )
  expect_type(conflicts_chr, "list")
  expect_length(conflicts_chr, 1)
  expect_equal(nrow(conflicts_chr[[1]]), nrow(conflicts_int[[1]]))
  expect_false(any(is.na(conflicts_int[[1]]$sample1)))
})

test_that("find_matching_rows, estimateLC, and methylation_score work across formats and edge cases", {
  # Probe ID matching for legacy, EPICv2, MSA, and MOUSE styles
  expect_equal(
    ewastools:::find_matching_rows(
      c("cg00000029", "cg00000108"),
      c("cg00000108_TC21", "cg00000029_BC11", "mu00000001_BC11", "rp00000002_TC11")
    ),
    c(2L, 1L)
  )

  k450 = testthat::test_path("../../data/5640269011_R01C01")
  stopifnot(file.exists(paste0(k450, "_Grn.idat")))

  beta = read_idats(k450, quiet = TRUE) |> correct_dye_bias() |> beta_values()

  # Leukocyte deconvolution (unconstrained and constrained)
  lc_uncon = estimateLC(beta, ref = "Reinius", constrained = FALSE)
  lc_con   = estimateLC(beta, ref = "Reinius", constrained = TRUE)
  expect_s3_class(lc_uncon, "data.table")
  expect_equal(nrow(lc_uncon), 1)
  expect_true(all(unlist(lc_uncon[1]) >= -1e-10))
  expect_equal(sum(unlist(lc_con[1])), 1, tolerance = 1e-6)

  # Methylation scores across all supported models
  for (mdl in c("cord_blood:maternal_smoking", "cord_blood:gestational_age", "placenta:gestational_age", "horvath_clock")) {
    sc = suppressWarnings(methylation_score(beta, mdl))
    expect_type(sc, "double")
    expect_length(sc, 1)
    expect_false(is.na(sc))
  }
  expect_null(methylation_score(beta, "nonexistent_model"))

  # Edge cases: zero matching probes or all-NA sample returns NA_real_ (not bare intercept)
  tiny_beta = matrix(c(0.5, NA_real_, 0.3, NA_real_), nrow = 2, ncol = 2, byrow = TRUE,
                     dimnames = list(c("cg00000029", "cg02256631"), c("ok", "all_na")))
  expect_warning(
    sc_partial <- methylation_score(tiny_beta, "cord_blood:maternal_smoking"),
    "not found"
  )
  expect_false(is.na(sc_partial["ok"]))
  expect_true(is.na(sc_partial["all_na"]))

  no_match_beta = tiny_beta[1, , drop = FALSE]
  expect_warning(
    sc_none <- methylation_score(no_match_beta, "horvath_clock"),
    "not found"
  )
  expect_equal(unname(sc_none), c(NA_real_, NA_real_))
})
