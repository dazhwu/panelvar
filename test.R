library(panelvar)
data(abdata)
p1 <- pvargmm(
  dependent_vars = c("n", "w", "k"),
  lags = 2,

  # exog_vars=c("k"),
  transformation = "fd",
  data = abdata,
  panel_identifier = c("id", "year"),
  steps = c("twostep"),
  system_instruments = TRUE,
  # max_instr_dependent_vars = 3,
  # max_instr_predet_vars = 3,
  # min_instr_dependent_vars = 1L,
  # min_instr_predet_vars = 1L,
  collapse = FALSE,
  progressbar = FALSE
)
summary(p1)

girf(p1, n.ahead = 8, ma_approx_steps = 8)
ex1_dahlberg_data_bs <- bootstrap_irf(p1,
  typeof_irf = c("GIRF"),
  n.ahead = 8,
  nof_Nstar_draws = 200,
  confidence.band = 0.95
)
print(ex1_dahlberg_data_bs)
