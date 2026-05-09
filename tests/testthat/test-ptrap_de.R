test_that("genes.filter returns only requested genes", {
  res <- tan_raw |>
    ptrap_de(treatment_name = "PACAP", genes.filter = c("Cdr1", "Pwwp2a"))
  expect_equal(nrow(res), 2L)
  expect_setequal(res$Gene, c("Cdr1", "Pwwp2a"))
})

test_that("genes.filter = NULL returns full result (no filtering)", {
  res_full <- tan_raw |> ptrap_de(treatment_name = "PACAP")
  res_nofilter <- tan_raw |>
    ptrap_de(treatment_name = "PACAP", genes.filter = NULL)
  expect_equal(nrow(res_full), nrow(res_nofilter))
})

test_that("genes.filter with a single gene returns one row", {
  res <- tan_raw |>
    ptrap_de(treatment_name = "PACAP", genes.filter = "Cdr1")
  expect_equal(nrow(res), 1L)
  expect_equal(res$Gene, "Cdr1")
})

test_that("genes.filter with non-existent gene returns 0 rows silently", {
  expect_no_warning(
    res <- tan_raw |>
      ptrap_de(
        treatment_name = "PACAP",
        genes.filter   = "ThisGeneDoesNotExist_xyz"
      )
  )
  expect_equal(nrow(res), 0L)
})

test_that("genes.filter does not affect kable.out path (returns knitr_kable)", {
  res <- tan_raw |>
    ptrap_de(
      treatment_name = "PACAP",
      genes.filter   = c("Cdr1", "Pwwp2a"),
      kable.out      = TRUE
    )
  expect_s3_class(res, "knitr_kable")
})
