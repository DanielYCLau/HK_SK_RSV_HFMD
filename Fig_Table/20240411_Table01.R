




#
#--- ESTIMATES OVER TIME ---
#






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

library("lubridate")





# version based on WIS table
sim.list = list()

for (i in 1:4) {
  # i = 1
  if (i == 1) {
    place   = "HK"
    virus   = "RSV"
    s.date  = seq(as.Date("2018-10-01"), as.Date("2023-10-01"), by = "6 months")
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
  version = "20240411"

  
  
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
  ( pmcmc.sim.path     = list.files("Fig_Table/sim", pmcmc.sim.filename, full.names = T) )
  load(pmcmc.sim.path)
  
  
  
  sim.i.list = lapply(
    c("s_all", "CCovid_p"), # , "Re"),
    function(par.k) {
      # par.k = "s"
      sim.i = subset(pmcmc.sim, time > 0)
      if (par.k == "CCovid_p") { sim.i[, par.k] = 1-sim.i[, par.k] } 
      sim.i = merge(x = date.seq, y = sim.i, by = "time", all.y = T)
      sim.i[, "week"] = week(sim.i[, "date"])
      
      sim.i.19 = subset(sim.i, year(date) == 2019)[, c(".id", "date", "week", par.k)]
      sim.i.20 = subset(sim.i, date >= as.Date("2020-01-01") & date < as.Date("2020-04-01")) [, c(".id", "date", "week", par.k)]
      sim.i.21 = subset(sim.i, date >= as.Date("2021-10-01") & date < as.Date("2022-04-01")) [, c(".id", "date", "week", par.k)]
      sim.i.22 = subset(sim.i, date >= as.Date("2022-04-01") & date < as.Date("2023-04-01")) [, c(".id", "date", "week", par.k)]
      sim.i.23 = subset(sim.i, date >= as.Date("2023-04-01") & date < as.Date("2024-03-01")) [, c(".id", "date", "week", par.k)]

      
      
      delta = lapply(
        sprintf("sim.i.%d", 20:23),
        function(y) {
          # y = "sim.i.20"
          data.y = get(y)
          colnames(sim.i.19)[4] = "par"
          colnames(data.y)[4] = "par"
          df.y = merge(x = sim.i.19, y = data.y, all.y = T, by = c(".id", "week") )
          df.y[, "delta"] = with(df.y, par.y - par.x)
          df.y = aggregate(delta ~ date.y, data = df.y, FUN = CI95.fx)
          df.y = do.call(data.frame, df.y)
          colnames(df.y) = c("time", "mean", "sd", "Q50", "LB", "UB", "min", "max")
          # df.y = rbind(df.y, NA)
          return(df.y)
        }
      )
      
      
      minmax = lapply(
        1:length(delta),
        function(p) {
          df = delta[[p]]
          if (par.k == "Re") {
            df[, "mean.95CrI"] = with(df, sprintf("%.1f (%.1f, %.1f)", mean, LB, UB))
          } else {
            df[, "mean.95CrI"] = with(df, sprintf("%.1f (%.1f, %.1f)", mean*100, LB*100, UB*100))
          }
          tmp = with(
            df, data.frame(
              "place"   = place,
              "virus"   = virus,
              "par"     = par.k,
              "ind"     = p,
              "period"  = paste0(format(range(time), "%b%y"), collapse = " - "),
              "min.t"   = time[which.min(mean)],
              "max.t"   = time[which.max(mean)],
              "min.est" = mean.95CrI[which.min(mean)],
              "max.est" = mean.95CrI[which.max(mean)]
            )
          )
          return(tmp)
        }
      )
      
      minmax = do.call(rbind, minmax)
      return(minmax)
    }
  )
  
  sim.i = do.call(rbind, sim.i.list)

  sim.list[[sprintf("%s_%s", place, virus)]] = sim.i
  
}





Tb1 = lapply(
  1:length(sim.list),
  function(k) {
    tmp.k = with(sim.list[[k]],  data.frame(par, ind, period, sprintf("%s - %s", min.est, max.est)))
    colnames(tmp.k)[4] = names(sim.list)[[k]]
    return(tmp.k)
  }
)

Tb1 = Reduce(function(x, y) merge(x, y, all = T), Tb1)
Tb1 = Tb1[with(Tb1, rev(order(par, -ind))),]







#
#--- OUTPUT ---
#





write.csv(
  Tb1,
  file      = sprintf("Fig_Table/%s_Table01.csv", format(Sys.Date(), "%Y%m%d")), 
  na        = "",
  row.names = F
)





#
#--- END ---
#




