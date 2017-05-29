test_that("catmap.funnel works", {
  expect_equal(data(catmapdata), )
  catmapobject <- catmap(catmapdata, 0.95, TRUE)
  expect_equal(catmap.funnel(catmapobject, TRUE), )
})
