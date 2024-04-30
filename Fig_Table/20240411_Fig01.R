




#
#--- FIG 01 ---
#





sim.filename = sprintf("%s_Fig01.pdf", format(Sys.Date(), "%Y%m%d"))
pdf(width = 8, height = 7, file = sprintf("Fig_Table/%s", sim.filename))

# windows(width = 10, height = 12)
# quartz(width = 10, height = 12)
par(mfcol = c(4,1), mar = c(0,7,0,5), oma = c(3,0,1,1), las = 1)



for (i in 1:4) {
  
  # i = 1
  if (i == 1) {
    place   = "HK"
    virus   = "RSV"
    y.lim   = 0.015
    leg.tx  = "A) HK-RSV"
  } else if (i == 2) {
    place   = "HK"
    virus   = "HFMD"
    y.lim   = 0.015
    leg.tx  = "B) HK-HFMD"
  } else if (i == 3) {
    place   = "SK"
    virus   = "RSV"
    y.lim   = 0.04
    leg.tx  = "C) SK-RSV"
  } else if (i == 4) {
    place   = "SK"
    virus   = "HFMD"
    y.lim   = 0.03
    leg.tx  = "D) SK-HFMD"
  }
  version = "20240411"
  
  
  setwd("program")
  source(sprintf("%s_SEIR_fx.R", version))
  proj.wd = getwd()
  
  data.formatting.fx(place, virus, scale = T)
  
  if (virus == "RSV") {
    data = RSV.data
  } else if (virus == "HFMD") {
    data = HFMD.data
  }
  
  
  
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
  pomp.data[, "dec.date"] = decimal_date(pomp.data[, "date"])
  
  
  
  
  # load sim
  pmcmc.sim.filename = sprintf("%s_%s_SEIR_%s_pmcmc_post_sim.Rdata", version, place, virus)
  ( pmcmc.sim.path     = list.files(sprintf("%s/../Fig_Table/sim", proj.wd), pmcmc.sim.filename, full.names = T) )
  load(pmcmc.sim.path)
  
  pmcmc.sim = merge(x = pmcmc.sim, y = pomp.data[, c("time", "pop", "all")], by = "time", all.x = T)
  pmcmc.sim = transform(pmcmc.sim, obs.p = obs / pop, obs.a = obs / all )
  
  
  
  
  # windows(width = 16, height = 12)
  plot(
    NULL, 
    xlim = c(2014, 2024),
    ylim = c(0, y.lim) * 1.25,
    xaxt = "n",
    yaxt = "n",
    ylab = ""
  )
  axis(1, at = 2014:2024, labels = F)
  if (i == 4) { 
    axis(1, at = c(2014:2023 + 0.5), labels = 2014:2023, tick = F)
  }
  abline(v = 2014:2024, col = alpha(1, 0.1), lty = 3, lwd = 1.5)
  if (i %in% 1:2) {
    axis(2, at = seq(0, 0.015, 0.005), labels = sprintf("%.3f", seq(0, 0.015, 0.005)) )
  } else if (i == 3) {
    axis(2, at = seq(0, 0.04, 0.01), labels = sprintf("%.3f", seq(0, 0.04, 0.01)) )
  } else if (i == 4) {
    axis(2, at = seq(0, 0.03, 0.01), labels = sprintf("%.3f", seq(0, 0.03, 0.01)) )
  }
  
  
  obs.i = aggregate(obs.p ~ time, data = pmcmc.sim, FUN = CI95.fx)
  obs.i = do.call(data.frame, obs.i)
  colnames(obs.i) = gsub("obs[.]p[.]", "", colnames(obs.i))
  obs.i = merge(x = obs.i, y = pomp.data, all.x = T)
  
  

  with(
    obs.i, {
      lines(mean ~ dec.date)
      polygon.fx(x = dec.date, LB = LB, UB = UB, col = alpha(1, 0.1))
      lines(obs/pop ~ dec.date, col = alpha(2, 0.75), lwd = 2)
    }
  )
  
  
  if (place == "SK") {
    lines(
      scaled.data ~ dec.date, 
      data = subset(data, dec.date > 2014 & WeekEnd <= stdt ), 
      col  = alpha(2, 0.5), 
      lwd  = 2,
      lty  = 2
    )
  }
  
  
  
  
  
  
  par(new = T)
  plot(
    NULL, 
    xlim = c(2014, 2024),
    ylim = c(-1, 1),
    xaxt = "n",
    yaxt = "n",
    ylab = ""
  )
  axis(4, at = seq(0, 1, 0.2), col = "#ff7f0e", col.axis = "#ff7f0e")
  
  lines(
    index/100 ~ dec.date, 
    data = subset(cov.mat, WeekEnd >= as.Date("2020-01-01") & WeekEnd <= as.Date("2023-02-28")), 
    lty = 2, 
    col = alpha("#1f77b4", 0.5)
  )

  covid_p.i = aggregate(CCovid_p ~ time, data = pmcmc.sim, FUN = CI95.fx)
  covid_p.i = do.call(data.frame, covid_p.i)
  colnames(covid_p.i) = gsub("CCovid_p[.]", "", colnames(covid_p.i))
  covid_p.i = merge(x = covid_p.i, y = pomp.data, all.x = T)
  
  lines(mean ~ dec.date, data = covid_p.i, col = alpha("#ff7f0e", 0.5))
  with(
    covid_p.i, 
    polygon.fx(x = dec.date, LB = LB, UB = UB, col = alpha("#ff7f0e", 0.1))
  )

  
  

    
  legend(
    "topleft",legend = leg.tx, bty = "n", cex = 1.5, text.font = 2
  )

  if (i == 3) {
    text(
      x = 2025.2, y = 1.5, labels = "Proportional Transmissibility", 
      col = "#ff7f0e", cex = 1.25, srt = 270, xpd = NA, font = 2
    )
    
  }
  
  setwd(sprintf("%s/../", proj.wd))
  
}

mtext("Virus activity (scaled proxy)", 2, las = 0, line = -2, outer = T, font = 2)





dev.off()










#
#--- END ---
#




