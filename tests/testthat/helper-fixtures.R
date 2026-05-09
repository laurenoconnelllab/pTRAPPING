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
