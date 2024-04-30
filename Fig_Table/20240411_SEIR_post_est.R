




#
#--- ESTIMATES ---
#





library(coda)
setwd("program")
proj.wd = getwd()



CI95.fx = function(x, round = F, dg = 2) {
  # summary statistics from samples 
  x = na.omit(x)
  sum.x = c(
    "mean"   = mean(x),
    "sd"     = sd(x),
    "median" = median(x),
    "LB"     = unname(quantile(x, 0.025)),
    "UB"     = unname(quantile(x, 0.975)),
    "min"    = min(x),
    "max"    = max(x)
  )
  if (round == T) {
    sum.x = round(sum.x, dg)
  }
  return(sum.x)
}

expit.fx = function(x) { exp(x)/(1 + exp(x)) }





# version based on WIS table
est.df = lapply(
  1:4, #c(1,3,2,4),
  function(i) {
    
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
    

    
    pmcmc.coda.filename = sprintf("%s_SEIR_%s_pmcmc_post2_coda_%02dc_%.fkiter.Rdata", place, virus, 60, iter)
    ( pmcmc.coda.filename = list.files(proj.wd, pmcmc.coda.filename, full.names = T) )
    load(pmcmc.coda.filename)
    

    
    pmcmc.coda = subset(pmcmc.coda, select = -c(.id, loglik, log.prior))
    
    est.i = lapply(
      colnames(pmcmc.coda),
      function(par.k) {
        coda.k = pmcmc.coda[, par.k]
        
        if (grepl("psi", par.k)) {
          new.k  = gsub("log_", "", par.k)
          coda.k = exp(coda.k)
        } else if (grepl("phi|logit", par.k)) {
          new.k  = gsub("logit_", "", par.k)
          coda.k = expit.fx(coda.k)
        } else if (i == 2 & par.k == "beta1") {
          # new.k  = "zeta1"
          new.k  = par.k # "beta1_v2"
          coda.k = exp(coda.k) + 0.2
        } else if (i == 4 & par.k == "beta2") {
          new.k  = par.k
          coda.k = -exp(coda.k)
        } else if (i == 1 & par.k == "b_c1") {
          new.k  = par.k # "beta1_v2"
          coda.k = -exp(coda.k)
        } else if (par.k == "b_c2") {
          new.k  = sprintf("%s_100", par.k)
          coda.k = exp(coda.k) * 100
        } else {
          new.k  = par.k
        }
        est.k = as.data.frame(t(CI95.fx(coda.k)))
        y = data.frame(
          "par"      = new.k,
          "mu.CI95"  = with(est.k, sprintf("%.2f (%.2f, %.2f)",   mean, LB, UB))
          # "Q50.CI95" = with(est.k, sprintf("%.2f (%.2f, %.2f)", median, LB, UB))
        )
        return(y)
      }
    )
    
    est.i = do.call(rbind, est.i)
    
    colnames(est.i)[2] = sprintf("%s.%s.%s", place, virus, colnames(est.i)[2])
    
    return(est.i)
  }
)





#
#--- OUTPUT ---
#





est.df = Reduce(function(x, y) merge(x, y, by = "par", all = T), est.df)

par.order = match(
  c("beta1", "beta2", "beta3", "beta4", "beta5",
    "phi", "psi_1",
    "seed", "b_c1", "b_c2_100", "b_c3", "psi_2"),
  est.df[, "par"]
)

est.df = est.df[par.order, ]        

write.csv(
  est.df, 
  file      = sprintf("../Fig_Table/%s_SEIR_est_table.csv", format(Sys.Date(), "%Y%m%d")), 
  na        = "-",
  row.names = F
)





#
#--- END ---
#




