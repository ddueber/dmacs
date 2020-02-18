test_that("lavaan_dmacs - Group - Continuous/Categorical", {

  cont_test_data <- read.csv("SampleData.csv")
  cat_test_data <- as.data.frame(lapply(cont_test_data, ordered))
  UniModel <- "F1 =~ F1_1 + F1_2 + F1_3 + F1_4 + F1_5 + F1_6 + F1_7 + F1_8 + F1_9 + F1_10"
  MultiModel <- "F1 =~ F1_1 + F1_2 + F1_3 + F1_4 + F1_5 + F1_6 + F1_7 + F1_8 + F1_9 + F1_10
                 F2 =~ F2_1 + F2_2 + F2_3 + F2_4 + F2_5 + F2_6 + F2_7 + F2_8 + F2_9 + F2_10 + F2_11 + F2_12 + F2_13 + F2_14 + F2_15 + F2_16 + F2_17 + F2_18 + F2_19 + F2_20
                 F3 =~ F3_1 + F3_2 + F3_3 + F3_4 + F3_5 + F3_6 + F3_7 + F3_8 + F3_9 + F3_10 + F3_11 + F3_12
                 F4 =~ F4_1 + F4_2 + F4_3 + F4_4 + F4_5"
  #fit0000 <- lavaan::cfa(model = MultiModel, data = cont_test_data, group = "Group")
  #fit0001 <- lavaan::cfa(model = MultiModel, data = cat_test_data, group = "Group")
  fit0010 <- lavaan::cfa(model = UniModel, data = cont_test_data, group = "Group")
  #fit0011 <- lavaan::cfa(model = UniModel, data = cat_test_data, group = "Group")
  fit0100 <- lavaan::cfa(model = MultiModel, data = cont_test_data[cont_test_data$Group < 3, ], group = "Group")
  #fit0101 <- lavaan::cfa(model = MultiModel, data = cat_test_data[cat_test_data$Group < 3, ], group = "Group")
  #fit0110 <- lavaan::cfa(model = UniModel, data = cont_test_data[cont_test_data$Group < 3, ], group = "Group")
  fit0111 <- lavaan::cfa(model = UniModel, data = cat_test_data[cat_test_data$Group < 3, ], group = "Group")

  load("lavsamples.RData")

  #expect_equal(lavaan_dmacs(fit0000, RefGroup = "2"), lavaan0000, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0001, RefGroup = "2"), lavaan0001, tolerance = .00001)
  expect_equal(lavaan_dmacs(fit0010, RefGroup = "2"), lavaan0010, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0011, RefGroup = "2"), lavaan0011, tolerance = .00001)
  expect_equal(lavaan_dmacs(fit0100, RefGroup = "2"), lavaan0100, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0101, RefGroup = "2"), lavaan0101, tolerance = .00001)
  #expect_equal(lavaan_dmacs(fit0110, RefGroup = "2"), lavaan0110, tolerance = .00001)
  expect_equal(lavaan_dmacs(fit0111, RefGroup = "2"), lavaan0111, tolerance = .00001)

})

test_that("lavaan_dmacs - Longitudinal - Continuous/Categorical", {

  cont_test_data <- read.csv("SampleDataLong.csv")
  cat_test_data <- as.data.frame(lapply(cont_test_data, ordered))
  LongModel <- "F1 =~ F1_1 + F1_2 + F1_3 + F1_4 + F1_5
                 F2 =~ F2_1 + F2_2 + F2_3 + F2_4 + F2_5
                 F3 =~ F3_1 + F3_2 + F3_3 + F3_4 + F3_5"
  fit1000 <- lavaan::cfa(model = LongModel, data = cont_test_data)
  fit1001 <- lavaan::cfa(model = LongModel, data = cat_test_data)

  load("lavsamplesLong.RData")

  expect_equal(lavaan_dmacs(fit1000, RefGroup = "F2", MEtype = "Long"), lavaan1000, tolerance = .00001)
  expect_equal(lavaan_dmacs(fit1001, RefGroup = "F2", MEtype = "long"), lavaan1001, tolerance = .00001)

})

