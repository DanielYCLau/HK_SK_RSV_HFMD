

# POST #


#
#--- PREPARATION ---
#





version = "20240411"
setwd("program")
source(sprintf("%s_SEIR_fx.R", version))
proj.wd = getwd()



place = "HK"
virus = "RSV"

data.formatting.fx(place, virus, scale = T)



if ( !dir.exists(pomp.dir) ) { dir.create(pomp.dir) }
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





#
#--- PMCMC ---
#





pmcmc.coda.filename = sprintf("%s_SEIR_%s_pmcmc_post_coda.*Rdata", place, virus)
load(list.files(proj.wd, pmcmc.coda.filename, full.names = T))



# win.fx(width = 12, height = 12); pairs(pmcmc.init[, -c(1,3)])

pmcmc.init = arrange(pmcmc.coda, -loglik)
pmcmc.init = pmcmc.init[!duplicated(pmcmc.init), ]
pmcmc.init = pmcmc.init[1:60, ]
pmcmc.init = subset(pmcmc.init, select = -c(.id, loglik, log.prior))

sigma      = cov(pmcmc.init)

par.proposal  = mvn.rw.adaptive(
  rw.var      = sigma,
  scale.start = 1000, 
  shape.start = 1000
)

# windows()
# pairs(pmcmc.init)



n.chain = 60
Npf  = 500
iter = 30000


# chain = 10; Npf = 500; iter =  1000; 2h
# chain = 60; Npf = 500; iter = 15000; 21h
# chain = 60; Npf = 500; iter = 30000; 40-50h

start.t = Sys.time()

pb = txtProgressBar(max = n.chain, style = 3)
progress = function(i) { setTxtProgressBar(pb, i) }
opts = list(progress = progress)

cl = makeCluster(n.chain, type = "SOCK")
registerDoSNOW(cl)

pmcmc.fit = foreach(
  i         = 1:n.chain,
  .combine  = c,
  .packages = "pomp",
  .inorder  = F,
  .options.snow = opts
) %dopar% {
  
  set.seed(i)
  
  pmcmc.init.i = pmcmc.init[i, ]
  
  pmcmc.fit.i = pmcmc(
    data     = pomp.model,
    Nmcmc    = iter,
    Np       = Npf,
    params   = pmcmc.init.i,
    proposal = par.proposal
  )
  
  # return(list(pmcmc.fit.i))
  trace.i = traces(pmcmc.fit.i)
  
  return(list(trace.i))
}

stopCluster(cl)

end.t = Sys.time()

cat("\n")
print(end.t - start.t)



pmcmc.fit.filename = file.fx(sprintf("%s_SEIR_%s_pmcmc_post2_fit_%02dc_%.fkiter.Rdata", place, virus, n.chain, iter/1e3))
save(pmcmc.fit, file = sprintf("%s/%s", proj.wd, pmcmc.fit.filename))
# load(list.files(proj.wd, "pmcmc_post2_fit.*Rdata", full.names = T))





#
#--- PMCMC result ---
#





n.chain = 60
iter    = 30000




# load data
pmcmc.fit.filename = sprintf("%s_SEIR_%s_pmcmc_post2_fit_%02dc_%.fkiter.Rdata", place, virus, n.chain, iter/1e3)
load(list.files(proj.wd, pmcmc.fit.filename, full.names = T))



#--- convergence diagnosis ---
pmcmc.trace = do.call(mcmc.list, pmcmc.fit)

pmcmc.trace = window(pmcmc.trace, start = 10000)
round( rejectionRate(pmcmc.trace), 2 )

# pmcmc.AC = autocorr.diag(pmcmc.trace, lags=c(1,10,20,50,100,200))
# apply(pmcmc.AC[, par.names], 2, round, digits = 3)

pmcmc.trace = window(pmcmc.trace, thin = 100)
# summary(pmcmc.trace)

# summary
( pmcmc.est = pmcmc.sum.fx(pmcmc.trace) )
pmcmc.est.filename = file.fx(sprintf("%s_SEIR_%s_pmcmc_post2_est.csv", place, virus))
write.csv(pmcmc.est, file = sprintf("%s/%s", proj.wd, pmcmc.est.filename), row.names = F)

# save coda 
pmcmc.coda = lapply(pmcmc.trace, as.data.frame)
pmcmc.coda = bind_rows(pmcmc.coda, .id = ".id")
pmcmc.coda = transform(pmcmc.coda, .id = as.numeric(as.character(.id)) )
rownames(pmcmc.coda) = NULL

pmcmc.coda.filename = file.fx(sprintf("%s_SEIR_%s_pmcmc_post2_coda_%02dc_%.fkiter.Rdata", place, virus, n.chain, iter/1e3))
save(pmcmc.coda, file = sprintf("%s/%s", proj.wd, pmcmc.coda.filename))





#
#--- END ---
#




