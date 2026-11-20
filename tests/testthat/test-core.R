test_that("out-of-sample projection reproduces the training scores", {
  set.seed(1); X <- matrix(rnorm(60 * 12), 60, 12)
  for (dc in c("pca", "cosine", "rbf", "poly")) {
    pr <- fit_projection(X, dc, d = 5); Z1 <- apply_projection(pr, X)
    expect_equal(dim(Z1), c(60, 5)); expect_true(max(abs(colMeans(Z1))) < 1e-8); expect_true(max(abs(apply(Z1, 2, sd) - 1)) < 1e-8)
    expect_equal(dim(apply_projection(pr, X[1:7, ])), c(7, 5))
  }
})

test_that("slides of a patient are averaged into one row", {
  sim <- simulate_cohort(20, D = 6, tiles_mean = 25, p_multi = 0.5, seed = 2)
  X <- pool_cohort(sim$tiles, sim$patient)
  expect_equal(nrow(X), length(unique(sim$patient))); expect_equal(rownames(X), unique(sim$patient))
  expect_length(pool_tiles(sim$tiles[[1]], "moment"), 12)
})

test_that("estimate, null and record behave", {
  sim <- simulate_hidden_context(n = 120, pi_A = 1, p = 8, seed = 3, mc = 0)
  f <- pathmi_fit(sim$Z, sim$y, nmc = 400, burn = 100, thin = 5)
  expect_true(f$estimate >= 0 && f$estimate <= f$entropy)
  nl <- pathmi_null(sim$Z, sim$y, perms = draw_permutations(120, 5), observed = f, nmc = 400, burn = 100, thin = 5)
  expect_length(nl$null, 5)
  hs <- half_sampling(sim$Z, sim$y, B = 4, nmc = 400, burn = 100, thin = 5)
  expect_equal(hs$se, sd(hs$draws) * sqrt(hs$m / (120 - hs$m)))
  eps <- reproducibility_bound(c(0, 0.001), c(0, -0.002))
  expect_equal(eps$bound, 0.002)
  rec <- recoverability_record(list(model = "simulated"), f, nl, hs, eps$bound)
  expect_s3_class(rec, "recoverability_record")
  skip_if_not_installed("jsonlite")
  js <- jsonlite::fromJSON(record_to_json(rec)); expect_equal(js$measurement$represented, f$estimate, tolerance = 1e-6)
})

test_that("the tile decomposition is exact", {
  sim <- simulate_cohort(80, D = 10, tiles_mean = 30, p_multi = 0.2, effect = 1, seed = 4)
  td <- tile_decomposition(sim$tiles, sim$patient, sim$y, d = 6, nmc = 400, burn = 100, thin = 5)
  expect_lt(td$max_identity_error, 1e-8)
  expect_true(all(td$tiles$surprise >= 0))
})
