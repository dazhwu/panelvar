library(panelvar)
data(abdata)
p1 <-pvargmm(
  dependent_vars = c("n","w","k"),
  lags = 2,
  
  #exog_vars=c("k"),
  transformation = "fd",
  data = abdata,
  panel_identifier = c("id", "year"),
  steps = c("twostep"),
  system_instruments = TRUE,
  #max_instr_dependent_vars = 3,
  #max_instr_predet_vars = 3,
  #min_instr_dependent_vars = 1L,
  #min_instr_predet_vars = 1L,
  collapse = FALSE,
  progressbar = FALSE
)
summary(p1)