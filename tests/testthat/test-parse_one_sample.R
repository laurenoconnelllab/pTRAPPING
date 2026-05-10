test_that(".parse_one_sample() handles T-R-F ordering", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("PACAP1IP", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "IP")

  r2 <- pTRAPPING:::.parse_one_sample("PACAP_1_IP", "IP", "INPUT")
  expect_equal(r2$treatment, "PACAP")
  expect_equal(r2$block, "1")
  expect_equal(r2$fraction, "IP")
})

test_that(".parse_one_sample() handles T-F-R ordering (current baseline)", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("PACAP_Input_1", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "INPUT")

  r2 <- pTRAPPING:::.parse_one_sample("PACAP_IP_1", "IP", "INPUT")
  expect_equal(r2$treatment, "PACAP")
  expect_equal(r2$block, "1")
  expect_equal(r2$fraction, "IP")
})

test_that(".parse_one_sample() handles R-T-F ordering", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("1_PACAP_IP", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "IP")

  r2 <- pTRAPPING:::.parse_one_sample("1_PACAP_INPUT", "IP", "INPUT")
  expect_equal(r2$treatment, "PACAP")
  expect_equal(r2$block, "1")
  expect_equal(r2$fraction, "INPUT")
})

test_that(".parse_one_sample() handles R-F-T ordering", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("1_IP_PACAP", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "IP")

  r2 <- pTRAPPING:::.parse_one_sample("1_INPUT_PACAP", "IP", "INPUT")
  expect_equal(r2$treatment, "PACAP")
  expect_equal(r2$block, "1")
  expect_equal(r2$fraction, "INPUT")
})

test_that(".parse_one_sample() handles F-T-R ordering", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("IP_PACAP_1", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "IP")

  r2 <- pTRAPPING:::.parse_one_sample("INPUT_PACAP_1", "IP", "INPUT")
  expect_equal(r2$treatment, "PACAP")
  expect_equal(r2$block, "1")
  expect_equal(r2$fraction, "INPUT")
})

test_that(".parse_one_sample() handles F-R-T ordering", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("IP_1_PACAP", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "IP")

  r2 <- pTRAPPING:::.parse_one_sample("INPUT_1_PACAP", "IP", "INPUT")
  expect_equal(r2$treatment, "PACAP")
  expect_equal(r2$block, "1")
  expect_equal(r2$fraction, "INPUT")
})

test_that(".parse_one_sample() 'in' is accepted as alias for INPUT", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("PACAP_in_1", "IP", "INPUT")
  expect_equal(r$fraction, "INPUT")

  r2 <- pTRAPPING:::.parse_one_sample("1_in_PACAP", "IP", "INPUT")
  expect_equal(r2$fraction, "INPUT")
})

test_that(".parse_one_sample() is case-insensitive for fraction keywords", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("PACAP_input_1", "IP", "INPUT")
  expect_equal(r$fraction, "INPUT")

  r2 <- pTRAPPING:::.parse_one_sample("PACAP_ip_1", "IP", "INPUT")
  expect_equal(r2$fraction, "IP")
})

test_that(".parse_one_sample() handles mixed separators", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("PACAP.IP-1", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "IP")
})

test_that(".parse_one_sample() INPUT does not corrupt IP token (no separators)", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("PACAP1INPUT", "IP", "INPUT")
  expect_equal(r$treatment, "PACAP")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "INPUT")
})

test_that(".parse_one_sample() errors when fraction is missing", {
  devtools::load_all(quiet = TRUE)

  expect_error(
    pTRAPPING:::.parse_one_sample("PACAP_1", "IP", "INPUT"),
    "Cannot identify IP/INPUT fraction"
  )
})

test_that(".parse_one_sample() errors when replicate is missing", {
  devtools::load_all(quiet = TRUE)

  expect_error(
    pTRAPPING:::.parse_one_sample("PACAP_IP", "IP", "INPUT"),
    "Cannot identify a replicate number"
  )
})

# --- region_names tests -------------------------------------------------------

test_that(".parse_one_sample() extracts region and cleans treatment (T-R-F-Region)", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("b1inputPOA", "IP", "INPUT", region_names = "POA")
  expect_equal(r$treatment, "b")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "INPUT")
  expect_equal(r$region, "POA")
})

test_that(".parse_one_sample() region matching is case-insensitive", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("a1ipSTR", "IP", "INPUT", region_names = "str")
  expect_equal(r$treatment, "a")
  expect_equal(r$block, "1")
  expect_equal(r$fraction, "IP")
  expect_equal(r$region, "str") # stored in user-supplied case
})

test_that(".parse_one_sample() region stays NA when region_names not matched", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("b1inputPOA", "IP", "INPUT", region_names = "STR")
  expect_true(is.na(r$region))
  expect_equal(r$treatment, "bPOA") # POA becomes part of treatment
})

test_that(".parse_one_sample() region is NA when region_names = NULL", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("b1inputPOA", "IP", "INPUT")
  expect_true(is.na(r$region))
  expect_equal(r$treatment, "bPOA")
})

test_that(".parse_one_sample() handles region with explicit separators", {
  devtools::load_all(quiet = TRUE)

  r <- pTRAPPING:::.parse_one_sample("b_1_input_POA", "IP", "INPUT", region_names = "POA")
  expect_equal(r$treatment, "b")
  expect_equal(r$region, "POA")
})

test_that(".build_sample_df_from_cols() adds region column when regions detected", {
  devtools::load_all(quiet = TRUE)

  cols <- c("b1inputPOA", "b1ipPOA", "a1inputPOA", "a1ipPOA")
  df <- pTRAPPING:::.build_sample_df_from_cols(cols, "IP", "INPUT", region_names = "POA")
  expect_true("region" %in% names(df))
  expect_true(all(df$region == "POA"))
  expect_setequal(df$treatment, c("a", "b"))
})

test_that(".build_sample_df_from_cols() omits region column without region_names", {
  devtools::load_all(quiet = TRUE)

  cols <- c("PACAP1IP", "PACAP1INPUT")
  df <- pTRAPPING:::.build_sample_df_from_cols(cols, "IP", "INPUT")
  expect_false("region" %in% names(df))
})

test_that(".build_sample_df_from_cols() handles mixed region/non-region columns", {
  devtools::load_all(quiet = TRUE)

  cols <- c("b1inputPOA", "b1ipPOA", "b1inputSTR", "b1ipSTR")
  df <- pTRAPPING:::.build_sample_df_from_cols(cols, "IP", "INPUT", region_names = "POA")
  expect_true("region" %in% names(df))
  expect_equal(sum(df$region == "POA", na.rm = TRUE), 2L)
  expect_equal(sum(is.na(df$region)), 2L)
})
