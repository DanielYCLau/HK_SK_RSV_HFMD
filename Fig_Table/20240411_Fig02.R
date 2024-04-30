




#
#--- FIG 02 ---
#





sim.filename = sprintf("%s_Fig02.pdf", format(Sys.Date(), "%Y%m%d"))
pdf(width = 10, height = 8, file = sprintf("Fig_Table/%s", sim.filename))

# windows(width = 10, height = 12)
# quartz(width = 12, height = 10)
par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(4,2,0,0), las = 1)


for (i in 1:4) {
  
  # i = 1
  if (i == 1) {
    place   = "HK"
    virus   = "RSV"
    leg.tx  = "A) HK-RSV"
  } else if (i == 2) {
    place   = "HK"
    virus   = "HFMD"
    leg.tx  = "B) HK-HFMD"
  } else if (i == 3) {
    place   = "SK"
    virus   = "RSV"
    leg.tx  = "C) SK-RSV"
  } else if (i == 4) {
    place   = "SK"
    virus   = "HFMD"
    leg.tx  = "D) SK-HFMD"
  }
  version = "20240411"
  
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
  pomp.data[, "dec.date"] = decimal_date(pomp.data[, "date"])
  
  
  
  
  # load sim
  pmcmc.sim.filename = sprintf("%s_%s_SEIR_%s_pmcmc_post_sim.Rdata", version, place, virus)
  ( pmcmc.sim.path     = list.files(sprintf("%s/../Fig_Table/sim", proj.wd), pmcmc.sim.filename, full.names = T) )
  load(pmcmc.sim.path)
  

  

    
  s.i = aggregate(s_all ~ time, data = pmcmc.sim, FUN = CI95.fx)
  s.i = do.call(data.frame, s.i)
  colnames(s.i) = gsub("s_all[.]", "", colnames(s.i))
  s.i = merge(x = s.i, y = pomp.data, all.x = T)
  
  
  plot(
    NULL,
    xlim = c(0, 1),
    ylim = c(0, ifelse(virus == "RSV", 0.3, 0.12)),
    xaxt = "n",
    yaxt = "n",
    ylab = ""
  )

  axis(1, at = (0:12)/12, labels = NA)
  abline(v = (0:12)/12, col = alpha(1, 0.2), lty = 3)
  axis(1, at = (0:11 + 0.5)/12, labels = substr(month.abb, 1, 1), tick = F)
  if (virus == "RSV") {
    axis(2, at = seq(0, 0.3, 0.05))
  } else {
    axis(2, at = seq(0, 0.12, 0.02))
  }
  
  tmp = split(s.i, year(s.i[, "date"]))
  for (j in names(tmp)) {
    within(
      tmp[[j]], {
        x.j = dec.date %% 1
        if (as.numeric(j) <= 2019) {
          lines(x = x.j, y = mean, lty = 2, lwd = 2, col = alpha(1, 0.2))
        } else if (as.numeric(j) == 2020) {
          lines(x = x.j, y = mean, lwd = 2, col = col.fx(1))
          polygon.fx(x.j, LB, UB, col = col.fx(1, 0.1))
        } else if (as.numeric(j) == 2021) {
          lines(x = x.j, y = mean, lwd = 2, col = col.fx(2))
          polygon.fx(x.j, LB, UB, col = col.fx(2, 0.1))
        } else if (as.numeric(j) == 2022) {
          lines(x = x.j, y = mean, lwd = 2, col = col.fx(3))
          polygon.fx(x.j, LB, UB, col = col.fx(3, 0.1))
        } else if (as.numeric(j) == 2023) {
          lines(x = x.j, y = mean, lwd = 2, col = col.fx(4))
          polygon.fx(x.j, LB, UB, col = col.fx(4, 0.1))
        }
      }
    )
  }
  
  text(x = -0.025, y = ifelse(virus == "RSV", 0.3, 0.12)*0.975, labels = leg.tx, pos = 4, cex = 1.25, font = 2)


  
  setwd(sprintf("%s/../", proj.wd))
  
}

mtext("Population Susceptibility", 2, las = 0, line = 0.25, outer = T, font = 2)
mtext("Calendar Month", 1, las = 0, line = 0.2, outer = T, font = 2)

legend(
  "bottomleft", 
  legend = c("Before 2019     ", "2020", "2021", "2022", "2023"), 
  inset  = c(-0.9, -0.375), 
  xpd    = NA, 
  lty    = c(2, 1, 1, 1, 1), 
  lwd    = 2,
  col    = c(alpha(1, 0.5), col.fx(1:4)),
  horiz  = T, 
  bty    = "n"
)



dev.off()










#
#--- END ---
#




