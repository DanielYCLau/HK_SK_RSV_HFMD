




#
#--- ESTIMATES OVER TIME ---
#





version = "20240411"
setwd("program")
source(sprintf("%s_SEIR_fx.R", version))
proj.wd = getwd()





sim.list = list()

for (i in 1:4) {
  
  if (i == 1) {
    place   = "HK"
    virus   = "RSV"
    s.date  = seq(as.Date("2018-10-01"), as.Date("2024-04-01"), by = "6 months")
  } else if (i == 2) {
    place   = "HK"
    virus   = "HFMD"
    s.date  = seq(as.Date("2018-01-01"), as.Date("2024-01-01"), by = "6 months")
  } else if (i == 3) {
    place   = "SK"
    virus   = "RSV"
    s.date  = seq(as.Date("2018-01-01"), as.Date("2024-01-01"), by = "year")
  } else if (i == 4) {
    place   = "SK"
    virus   = "HFMD"
    s.date  = seq(as.Date("2018-01-01"), as.Date("2024-01-01"), by = "year")
  }

  
  
  if (place == "HK") {
    stdt = "2014-01-01"
  } else if (virus == "RSV") {
    stdt = "2016-10-01"
  } else if (virus == "HFMD") {
    stdt = "2016-01-01"
  }
  eddt = "2023-12-31"
  
  date.seq = seq(as.Date(stdt), as.Date(eddt), by = "day")
  date.seq = date.seq[weekdays(date.seq) == "Saturday"]
  date.seq = data.frame("date" = date.seq, "time" = 1:length(date.seq))
  

  
  
  
  # load sim
  pmcmc.sim.filename = sprintf("%s_%s_SEIR_%s_pmcmc_post_sim.Rdata", version, place, virus)
  ( pmcmc.sim.path     = list.files(sprintf("%s/../Fig_Table/sim", proj.wd), pmcmc.sim.filename, full.names = T) )
  load(pmcmc.sim.path)

  
  
  sim.i.list = lapply(
    c("s_all", "CCovid_p", "RR0", "RRe"),
    function(par.k) {
      
      sim.i = subset(pmcmc.sim, time > 0)
      if (par.k == "Covid_p") { sim.i[, par.k] = 1-sim.i[, par.k] } 
      sim.i = aggregate(get(par.k) ~ time, data = sim.i, FUN = CI95.fx)
      sim.i = do.call(data.frame, sim.i)
      colnames(sim.i) = c("time", "mean", "sd", "Q50", "LB", "UB", "min", "max")
      
      sim.i = merge(x = date.seq, y = sim.i, by = "time", all.y = T)
      sim.i = transform(sim.i, season = cut(date, s.date, right = F) )

      sim.i.y = split(sim.i, sim.i[, "season"])
      
      sim.i.y = lapply(
        names(sim.i.y),
        function(s) {
          tmp   = sim.i.y[[s]]
          if (nrow(tmp) > 0) {
            min.x = with(tmp, tmp[which.min(mean), ])
            max.x = with(tmp, tmp[which.max(mean), ])

            
            min.t = min.x[, "date"]
            max.t = max.x[, "date"]
            if (par.k %in% c("RR0", "RRe")) {
              min.k = with(min.x, sprintf("%.2f (%.2f, %.2f)", mean, LB, UB))
              max.k = with(max.x, sprintf("%.2f (%.2f, %.2f)", mean, LB, UB))
            } else {
              min.k = with(min.x, sprintf("%1.f%% (%1.f%%, %1.f%%)", mean*100, LB*100, UB*100))
              max.k = with(max.x, sprintf("%1.f%% (%1.f%%, %1.f%%)", mean*100, LB*100, UB*100))
            }
            
            
            y = data.frame(
              "place"  = place,
              "virus"  = virus,
              "par"    = par.k,
              "season" = format(as.Date(s), "%Y%b"),
              "min.t"  = min.t, "min"    = min.k,
              "max.t"  = max.t, "max"    = max.k
            )
            
            return(y)
            
          }
        }
      )
      
      sim.i.k = do.call(rbind, sim.i.y)
      sim.i.k[nrow(sim.i.k)+1, ] = NA
      return(sim.i.k)
    }
  )
  
  sim.i = do.call(rbind, sim.i.list)
  # sim.i[nrow(sim.i)+1, ] = NA
  
  sim.list[[i]] = sim.i
  
}





#
#--- OUTPUT ---
#





sim.est = do.call(rbind, sim.list)
write.csv(
  sim.est, 
  file      = sprintf("../Fig_Table/%s_SEIR_sim_est.csv", format(Sys.Date(), "%Y%m%d")), 
  na        = "",
  row.names = F
)





#
#--- END ---
#




