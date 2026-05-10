# Count matrix with brain-region-embedded column names (POA and STR).
# Treatment "a" and "b", 2 replicates each, both regions.
region_raw <- data.frame(
  Gene       = c("Adcyap1", "Bdnf", "Ucn3"),
  b1inputPOA = c(71L,  2L,  98L),
  b1ipPOA    = c(66L, 73L, 622L),
  b2inputPOA = c(60L, 25L,  28L),
  b2ipPOA    = c(440L, 325L, 930L),
  a1inputPOA = c(11L, 13L,  21L),
  a1ipPOA    = c(155L, 38L, 294L),
  a2inputPOA = c(100L, 74L,  36L),
  a2ipPOA    = c(277L, 54L, 133L),
  b1inputSTR = c(71L,  2L,  82L),
  b1ipSTR    = c(66L, 73L, 840L),
  b2inputSTR = c(60L, 25L,  51L),
  b2ipSTR    = c(440L, 325L, 1177L),
  a1inputSTR = c(11L, 13L,  93L),
  a1ipSTR    = c(155L, 38L, 765L),
  a2inputSTR = c(100L, 74L,  99L),
  a2ipSTR    = c(277L, 54L, 499L)
)

# Minimal count matrix used across ptrap_de tests.
# Columns are auto-parseable: treatment=PACAP, replicate=1/2/3, fraction=IP/INPUT.
tan_raw <- data.frame(
  Gene        = c("Adcyap1", "Bdnf", "Cdr1", "Gng8", "Nxph4", "Pwwp2a", "Ucn3"),
  PACAP1INPUT = c(200L, 150L,  80L, 300L, 250L, 120L, 180L),
  PACAP1IP    = c(800L, 600L, 320L, 1200L, 1000L, 480L, 720L),
  PACAP2INPUT = c(180L, 140L,  75L, 280L, 240L, 110L, 170L),
  PACAP2IP    = c(720L, 560L, 300L, 1120L,  960L, 440L, 680L),
  PACAP3INPUT = c(210L, 160L,  85L, 320L, 260L, 130L, 190L),
  PACAP3IP    = c(840L, 640L, 340L, 1280L, 1040L, 520L, 760L)
)
