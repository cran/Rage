test_that("qsd_converge works correctly", {
  x1 <- qsd_converge(mat_u)
  expect_length(x1, 1L)
  expect_gt(x1, 0)

  # matU all zero
  x2 <- suppressWarnings(qsd_converge(mat_u_zero))
  expect_identical(x2, 1L)

  # multiple starting stages
  x3 <- qsd_converge(mat_u, start = c(0.8, 0.2, 0, 0))
  expect_length(x3, 1L)
  expect_gt(x3, 0)

  # named life stages
  x4 <- qsd_converge(mat_u_named, start = "sm")
  expect_identical(x1, x4)
})

test_that("qsd_converge warns and fails gracefully", {
  expect_error(qsd_converge(mat_u_na))
  expect_error(qsd_converge(mat_u, start = 10))
  expect_error(qsd_converge(mat_u, start = "stage name"))
  expect_error(qsd_converge(mat_u_named, start = "invalid name"))
})

test_that("qsd_converge works w/ non-ergodic matrix", {
  mat_no_ergo <- matrix(
    c(
      0, 2.2, 5.8, 0, 0,
      0.8, 0, 0, 0, 0,
      0, 0.89, 0, 0, 0,
      0, 0, 0.93, 0, 0,
      0, 0, 0, 0.75, 0.5
    ),
    nrow = 5, ncol = 5, byrow = TRUE
  )

  # fast convergence
  t_qsd <- qsd_converge(mat_no_ergo)

  expect_length(t_qsd, 1L)

  # multi-state w/ non-ergodic also works.
  f_qsd <- qsd_converge(mat_no_ergo,
    start = c(1, 2, 3, 2, 0)
  )

  expect_length(f_qsd, 1L)
  # This depends on default value of conv in qsd_converge
})
