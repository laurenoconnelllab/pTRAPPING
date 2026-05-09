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
