
test_that("item dmacs computation is accurate", {
  expect_equal(item_dmacs(0.76, 0.65, 0.74, 1.28, 0.21, 1.76, 1.85), 0.3899706, tolerance = .00001)
  expect_equal(item_dmacs(0.76, 0.65, 0.74, 1.28, 0.21, 1.76, 1.85, TRUE), 0.08332242, tolerance = .00001)
  expect_equal(item_dmacs(0.54, c(0, 0.65), 0.89, c(-1.25, 1.28), 0.21, 1.76, 2.11), 0.06158231, tolerance = .00001)
  expect_equal(item_dmacs(-1.21, 0.71, -0.53, 1.11, 0.21, 1.76, 2.11), 0.5747149, tolerance = .00001)
  expect_error(item_dmacs(-1.21, 0.71, -0.53, c(-1.25, 1.28), 0.21, 1.76, 2.11), "Item must have same number of thresholds in both reference and focal group")
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

test_that("dmacs_summary is not broken", {

  LambdaList <- list(Group1 <- matrix(c(1.00, 0.74,  1.14, 0.92), ncol = 1),
                     Group2 <- matrix(c(1.00, 0.76,  1.31, 0.98), ncol = 1))
  LambdaList2 <- list(Group1 <- matrix(c(1.00, 0.74, 0, 0, 0, 0,  1.14, 0.92), ncol = 2),
                     Group2 <- matrix(c(1.00, 0.76, 0, 0, 0, 0,  1.31, 0.98), ncol = 2))
  ThreshList <- list(Group1 <- c(0.00, 1.28, -0.82, 0.44),
                     Group2 <- c(0.00, 0.65, -0.77, 0.47))
  ThreshList2 <- list(Group1 <- list(0.00, 1.28, -0.82, 0.44),
                     Group2 <- list(0.00, 0.65, -0.77, 0.47))
  ThreshList3 <- list(Group1 <- list(c(-1, 0.00, 1.5), c(-0.25, 1.28), c(-1.61, -0.82, 0, 1.11), 0.44),
                      Group2 <- list(c(-1.28, 0.23, 1.25), c(-0.25, 1.66), c(-1.61, -1.18, -0.25, 0.87), 0.25))
  ThreshList4 <- list(Group1 <- list(c(-1, 0.00, 1.5), c(-0.25, 1.28), c(-1.61, -0.82, 0, 1.11), 0.44),
                      Group2 <- list(c(-1.28, 0.23, 1.25), c(-0.25, 1.66), c(-1.61, -1.18, -0.25), 0.25))
  MeanList   <- list(Group1 <- 0.21,
                     Group2 <- 0.19)
  VarList    <- list(Group1 <- 1.76,
                     Group2 <- 1.34)
  SDList     <- list(Group1 <- c(2.12, 1.85,  1.12, 3.61),
                     Group2 <- c(2.38, 1.44,  1.67, 2.15))
  Groups <- c("Group1", "Group2")
  RefGroup <- "Group2"

  expect_equal(dmacs_summary(LambdaList, ThreshList, MeanList, VarList, SDList, Groups, RefGroup), list(DMACS = c(0.00000000, 0.38997060, 0.24811349, 0.02880497), ItemDeltaMean = c(0.00000000, 0.83021704, -0.11369383, -0.05651525), MeanDiff = 0.660008, VarDiff = -1.782), tolerance = .00001)
  expect_known_output(dmacs_summary(LambdaList2, ThreshList, MeanList, VarList, SDList, Groups, RefGroup), "d_sum_multi.rds", update = FALSE)
  expect_error(dmacs_summary(LambdaList, ThreshList, MeanList, VarList, SDList, Groups, RefGroup, categorical = TRUE),"Thresholds must be in a list indexed by item. The thresholds for each item should be a vector")
  expect_equal(dmacs_summary(LambdaList, ThreshList2, MeanList, VarList, SDList, Groups, RefGroup, categorical = TRUE), list(DMACS = c(0.000000000, 0.083322420, 0.022104064, 0.004459561), ItemDeltaMean = c(0.000000000, -0.160378048,  0.003161234,  0.011124317), MeanDiff = -0.1460925), tolerance = .00001)
  expect_equal(dmacs_summary(LambdaList, ThreshList3, MeanList, VarList, SDList, Groups, RefGroup), list(DMACS = c(0.02577838, 0.04575298, 0.22451775, 0.01621731), ItemDeltaMean = c(-0.05325132,  0.08741389, -0.26369274, -0.05811619), MeanDiff = -0.2876464), tolerance = .00001)
  expect_error(dmacs_summary(LambdaList, ThreshList4, MeanList, VarList, SDList, Groups, RefGroup), "Item must have same number of thresholds in both reference and focal group")
  expect_equal(dmacs_summary(LambdaList, ThreshList, MeanList, VarList, SDList), list(DMACS = c(0.00000000, 0.46819121, 0.13742292, 0.04046284), ItemDeltaMean = c(0.00000000, -0.72487849, 0.09526908, 0.04792394), MeanDiff = -0.5816855, VarDiff = 1.273), tolerance = .00001)
})

## dmacs_summary_single is tested by dmacs_summary

## lavaan_dmacs could use a more comprehensive suite of tests [all 8 versions!!]
test_that("lavaan_dmacs behaves appropriately", {

  HS.model <- '  visual  =~ x1 + x2 + x3
                 textual =~ x4 + x5 + x6
                 speed   =~ x7 + x8 + x9 '
  fit <- lavaan::cfa(HS.model,
                     data = lavaan::HolzingerSwineford1939,
                     group = "school")

  expect_known_output(lavaan_dmacs(fit, RefGroup = "Pasteur"), "HS_lav_dmacs.rds", update = FALSE)


  cont_test_data <- read.csv("SampleData.csv")
  cat_test_data <- as.data.frame(lapply(cont_test_data, ordered))
  UniModel <- "F1 =~ F1_1 + F1_2 + F1_3 + F1_4 + F1_5 + F1_6 + F1_7 + F1_8 + F1_9 + F1_10"
  MultiModel <- "F1 =~ F1_1 + F1_2 + F1_3 + F1_4 + F1_5 + F1_6 + F1_7 + F1_8 + F1_9 + F1_10
                 F2 =~ F2_1 + F2_2 + F2_3 + F2_4 + F2_5 + F2_6 + F2_7 + F2_8 + F2_9 + F2_10 + F2_11 + F2_12 + F2_13 + F2_14 + F2_15 + F2_16 + F2_17 + F2_18 + F2_19 + F2_20
                 F3 =~ F3_1 + F3_2 + F3_3 + F3_4 + F3_5 + F3_6 + F3_7 + F3_8 + F3_9 + F3_10 + F3_11 + F3_12
                 F4 =~ F4_1 + F4_2 + F4_3 + F4_4 + F4_5"
  fit0000 <- lavaan::cfa(model = MultiModel, data = cont_test_data, group = "Group")
  #fit0001 <- lavaan::cfa(model = MultiModel, data = cat_test_data, group = "Group")
  fit0010 <- lavaan::cfa(model = UniModel, data = cont_test_data, group = "Group")
  #fit0011 <- lavaan::cfa(model = UniModel, data = cat_test_data, group = "Group")
  fit0100 <- lavaan::cfa(model = MultiModel, data = cont_test_data[cont_test_data$Group < 3, ], group = "Group")
  #fit0101 <- lavaan::cfa(model = MultiModel, data = cat_test_data[cat_test_data$Group < 3, ], group = "Group")
  #fit0110 <- lavaan::cfa(model = UniModel, data = cont_test_data[cont_test_data$Group < 3, ], group = "Group")
  fit0111 <- lavaan::cfa(model = UniModel, data = cat_test_data[cat_test_data$Group < 3, ], group = "Group")

  load("lavsamples.RData")

  expect_equal(lavaan_dmacs(fit0000, RefGroup = "2"), lavaan0000, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0001, RefGroup = "2"), lavaan0001, tolerance = .00001)
  expect_equal(lavaan_dmacs(fit0010, RefGroup = "2"), lavaan0010, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0011, RefGroup = "2"), lavaan0011, tolerance = .00001)
  expect_equal(lavaan_dmacs(fit0100, RefGroup = "2"), lavaan0100, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0101, RefGroup = "2"), lavaan0101, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0110, RefGroup = "2"), lavaan0110, tolerance = .00001)
  expect_equal(lavaan_dmacs(fit0111, RefGroup = "2"), lavaan0111, tolerance = .00001)

})


test_that("mplus_dmacs behaves appropriately", {

  load("mplussamples.RData")

  #expect_equal(mplus_dmacs("mplus1000.out", RefGroup = "GROUP2"), mplus1000, tolerance = .00001)
  expect_equal(mplus_dmacs("mplus1001.out", RefGroup = "GROUP2"), mplus1001, tolerance = .00001)
  ## MplusAutomation cannot read sampstat from mplus1010.out for some completely unknown reason
  ##expect_equal(mplus_dmacs("mplus1010.out", RefGroup = "GROUP2"), mplus1010)
  expect_equal(mplus_dmacs("mplus1011.out", RefGroup = "GROUP2"), mplus1011, tolerance = .00001)
  #expect_equal(mplus_dmacs("mplus1100.out", RefGroup = "GROUP2"), mplus1100, tolerance = .00001)
  expect_equal(mplus_dmacs("mplus1101.out", RefGroup = "GROUP2"), mplus1101, tolerance = .00001)
  expect_equal(mplus_dmacs("mplus1110.out", RefGroup = "GROUP2"), mplus1110, tolerance = .00001)
  #expect_equal(mplus_dmacs("mplus1111.out", RefGroup = "GROUP2"), mplus1111, tolerance = .00001)


})
