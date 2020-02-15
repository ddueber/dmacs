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
