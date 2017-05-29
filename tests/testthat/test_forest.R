test_that("catmap.cumulative works", {
  expect_equal(data(catmapdata), )
  catmapobject <- catmap(catmapdata, 0.95, TRUE)
  expect_equal(catmap.cumulative(catmapobject, TRUE, TRUE, TRUE, TRUE), )
})

test_that("catmap.forest works", {
  expect_equal(data(catmapdata), )
  catmapobject <- catmap(catmapdata, 0.95, TRUE)
  expect_equal(catmap.forest(catmapobject, TRUE, TRUE), )
})

test_that("catmap.sense works", {
  expect_equal(data(catmapdata), )
  catmapobject <- catmap(catmapdata, 0.95, TRUE)
  expect_equal(catmap.sense(catmapobject, TRUE, TRUE, TRUE, TRUE), )
})
