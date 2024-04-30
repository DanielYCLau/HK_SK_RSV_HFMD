




#
#--- SIMULATION ---
#





start.t = Sys.time()

for (i in 1:4) {
  
  # i = 1
  if (i == 1) {
    place   = "HK"
    virus   = "RSV"
  } else if (i == 2) {
    place   = "HK"
    virus   = "HFMD"
  } else if (i == 3) {
    place   = "SK"
    virus   = "RSV"
  } else if (i == 4) {
    place   = "SK"
    virus   = "HFMD"
  }
  version = "20240411"
  iter    = 30
  
  setwd("program")
  source(sprintf("%s_SEIR_fx.R", version))
  proj.wd = getwd()
  
  data.formatting.fx(place, virus, scale = T)
  
  
  
  setwd(pomp.dir)
  pomp.dir = getwd()
  
  
  
  if (place == "HK") {
    stdt = "2014-01-01"
  } else if (virus == "RSV") {
    stdt = "2016-10-01"
  } else if (virus == "HFMD") {
    stdt = "2016-01-01"
  }
  eddt = "2023-12-31"
  pomp.model = model.spec.fx(version, place, virus, stdt, eddt)
  
  
  
  
  # load sim
  pmcmc.sim.filename = sprintf("%s_%s_SEIR_%s_pmcmc_post_sim.Rdata", version, place, virus)
  pmcmc.sim.path     = list.files(proj.wd, pmcmc.sim.filename, full.names = T)
  
  pmcmc.coda.filename = sprintf("%s_SEIR_%s_pmcmc_post2_coda_%02dc_%.fkiter.Rdata", place, virus, 60, iter)
  ( pmcmc.coda.filename = list.files(proj.wd, pmcmc.coda.filename, full.names = T) )
  load(pmcmc.coda.filename)
    
  pmcmc.sim = pmcmc.sim.fx(pmcmc.coda)
  pmcmc.sim = subset(pmcmc.sim, time > 0)
  save(pmcmc.sim, file = sprintf("%s/../Fig_Table/sim/%s", proj.wd, pmcmc.sim.filename))

  setwd(sprintf("%s/../", proj.wd))
  
}

end.t = Sys.time()

cat("\n")
print(end.t - start.t)





#
#--- END ---
#




