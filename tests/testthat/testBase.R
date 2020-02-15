
test_that("item dmacs computation is accurate", {
  expect_equal(item_dmacs(0.76, 0.74, 0.65, 1.28, 0.21, 1.76, 1.85), 0.3899707, tolerance = .00001)
  expect_equal(item_dmacs(0.76, 0.74, 0, 0, 0.21, 1.76, 1.85, 0.65, 1.28, 1, 1, categorical = TRUE), 0.1087381, tolerance = .00001)
  expect_equal(item_dmacs(0.54, 0.74, 1.21, 1.76, 2.11, 1.65, 1.62, c(-1.25, 1.28), c(-1, 0.65), 1, 1), 0.1417415, tolerance = .00001)
  expect_error(item_dmacs(0.54, 0.74, 1.21, 1.76, 2.11, 1.65, 1.62, c(-1.25, 1.28), 0.65, 1, 1), "Item must have same number of thresholds in both reference and focal group")
})

test_that("delta item mean computation is accurate", {
  expect_equal(delta_mean_item(0.76, 0.65, 0.74, 1.28, 0.21, 1.76), 0.830217, tolerance = .00001)
  expect_equal(delta_mean_item(0.76, 0.65, 0.74, 1.28, 0.21, 1.76, TRUE), -0.160378, tolerance = .00001)
  expect_equal(delta_mean_item(0.54, c(0, 0.65), 0.89, c(-1.25, 1.28), 0.21, 1.76), 0.1438971, tolerance = .00001)
  expect_equal(delta_mean_item(-1.21, 0.71, -0.53, 1.11, 0.21, 1.76), 0.7201052, tolerance = .00001)
  expect_error(delta_mean_item(-1.21, 0.71, -0.53, c(-1.25, 1.28), 0.21, 1.76), "Item must have same number of thresholds in both reference and focal group")
})

test_that("delta variance computation is accurate", {
  expect_equal(delta_var(c(1.00, 0.74,  1.14, 0.92), c(1.00, 0.76,  1.31, 0.98), 1.76), 1.672, tolerance = .00001)
  expect_null(suppressWarnings(delta_var(c(1.00, 0.74,  1.14, 0.92), c(1.00, 0.76,  1.31, 0.98), 1.76, categorical = TRUE)))
  expect_warning(delta_var(c(1.00, 0.74,  1.14, 0.92), c(1.00, 0.76,  1.31, 0.98), 1.76, categorical = TRUE), "Delta variance can only be computed for linear models, not for categorical ones")
})


test_that("expected value computation is accurate", {
  expect_equal(expected_value(1.21, -.30, .22, FALSE), -0.0338, tolerance = .00001)
  expect_equal(expected_value(1.21, c(-2.12, -.30, 1.80), .22, TRUE), 1.761025, tolerance = .00001)
  expect_equal(expected_value(1.21, c(-2.12, -.30, 1.80), .22, FALSE), 1.761025, tolerance = .00001)
})
