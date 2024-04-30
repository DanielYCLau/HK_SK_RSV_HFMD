




#
#--- PLOT SIM UNDER PROP. TRANSMISSIBILITY ---
#




obs.list = list()





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
  
  sen.model  = scan(sprintf("%s_SEIR_fx.R", version), what = character(), sep = "\n")
  sen.model  = sen.model[grep("SEIR.model = Csnippet", sen.model):(grep("rInit = Csnippet", sen.model)-1)]
  sen.model  = gsub("[(]b_c1_ind[)][?] -exp[(]b_c1[)] [:]", "", sen.model )
  SEIR.model = eval( parse(text = sen.model) )
  
  
  
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
  
  
  
  
  
  pmcmc.coda.filename = sprintf("%s_SEIR_%s_pmcmc_post2_coda_%02dc_%.fkiter.Rdata", place, virus, 60, iter)
  ( pmcmc.coda.filename = list.files(proj.wd, pmcmc.coda.filename, full.names = T) )
  load(pmcmc.coda.filename)
  

  
  
  
  # 10c; seq(-3, 1, 0.1); 50mins
  cl = makeCluster(min(n.cores, 11), type = "SOCK")
  registerDoSNOW(cl)

  if (place == "HK" & virus == "RSV") {
    sen.i = logit(seq(0.0, 1, length.out = 41))
  } else if (place == "SK" & virus == "RSV") {
    sen.i = logit(seq(0.5, 1, length.out = 41))
  } else {
    sen.i = logit(seq(0.7, 1, length.out = 41))
  }

  obs.i = foreach(
    sen       = sen.i,
    .packages = "pomp",
    .combine  = "rbind"
  ) %dopar% {
    coda.j = transform(pmcmc.coda, b_c1 = sen)

    sim.j = pmcmc.sim.fx(coda.j)#[1:100,])
    sim.j = subset(sim.j, time > 0)
    obs.j = aggregate(obs ~ time, data = sim.j, FUN = mean)
    obs.j[, "j"] = expit(sen)
    return(obs.j)
  }
  
  stopCluster(cl)
  


  
  
  # obs.i = do.call(rbind, obs.i)
  obs.i = merge(x = obs.i, y = pomp.data[, c("date", "time")], all.x = T)
  obs.i[, "dec.date"] = decimal_date(obs.i[, "date"])
  obs.i = subset(obs.i, dec.date >= 2019)
  obs.list[[sprintf("%s_%s", place, virus)]] = obs.i
  

  
  
  setwd(sprintf("%s/../", proj.wd))
  
}

end.t = Sys.time()
end.t - start.t





save(
  obs.list,
  file = sprintf("Fig_Table/sen_sim/%s", file.fx("SEIR_sen_sim.Rdata"))
)





#
#--- RESULT ---
#





library("lubridate")
library("lattice")
library("RColorBrewer")
library("scales")
library("pomp")

load(list.files("Fig_Table/sen_sim/", "sen_sim.Rdata", full.names = T))

for (i in 1:length(obs.list)) {
  obs.i = obs.list[[i]]
  max.i = max( subset(obs.i, dec.date < 2020)[, "obs"] )
  obs.i[, "obs"] = obs.i[, "obs"] / max.i
  obs.list[[i]]  = obs.i
}







# heat.col = colorRampPalette(rev(c("darkred","darkred","darkred", "red", "orange","yellow", "lightyellow")))
heat.col = colorRampPalette(rev(c("#510000", "#510000","darkred","darkred", "red", "orange", "yellow", "lightyellow", "white")))

main.par = list(
  par.main.text = list(
    font = 2,
    just = "top", 
    y    = grid::unit(-0.5, "mm")
    # x = grid::unit(5, "mm")
  )
)

for (i in 1:length(obs.list)) {
  
  assign(
    sprintf("%s_sim", names(obs.list)[i]),
    levelplot(
      obs ~ dec.date * j, 
      data = obs.list[[i]],  
      main = c("A) HK-RSV","B) HK-HFMD", "C) SK-RSV", "D) SK-HFMD")[i],
      par.settings = main.par,
      xlim = c(2019, decimal_date(as.Date("2023-02-28"))),
      xlab = "",
      ylab = expression(expit(r[0])),
      col.regions = heat.col,
      at          = seq(0, 2, length.out = 101),
      colorkey    = list(
        labels = list(
          at     = seq(0, 2, 0.5),
          labels = c(sprintf("   %.1f", seq(0, 1.5, 0.5)), expression("">="2.0"))
        )
      ),
      panel = function(...) {
        panel.fill(col = "darkred") 
        panel.levelplot(...)
        panel.abline(v = 2019:2023, lty = 3, col = alpha(1, 0.2))
        panel.abline(h = expit(c(-0.01, 1.19, 0.26, 1.82))[i], lty = 2, col = alpha(4, 0.5))
        # panel.abline(h = 0, lty = 2, col = alpha(4, 0.5))
      },
      xscale.components = function(...) {
        ans = xscale.components.default(...)
        ans$top = ans$bottom
        ans$top$ticks$at = 2019:2023
        ans$bottom$ticks$at = 2019:2023
        ans
      },
      scales = list(
        x = list(at = c(2019:2022 + 0.5, 2023.1), labels = 2019:2023),
        y = list(at = seq(0, 1, 0.1))
      )
    )
  )
  
}



pdf(
  height = 8, width = 12,
  file = sprintf("Fig_Table/%s_FigS03.pdf", format(Sys.Date(), "%Y%m%d"))
)

# windows(width = 12, height = 8)
# quartz(width = 12, height = 8)
lw = list(
  left.padding   = list(x =  0.1, units = "inches"),
  right.padding  = list(x =  0.1, units = "inches")
)
lh = list(
  top.padding    = list(x =  0.1, units = "inches"),
  bottom.padding = list(x = -0.2, units = "inches")
)

lattice.options(layout.widths = lw, layout.heights = lh)

i = 1; print(HK_RSV_sim,  split = c(1,1,2,2))
i = 2; print(HK_HFMD_sim, split = c(2,1,2,2), newpage = F)
i = 3; print(SK_RSV_sim,  split = c(1,2,2,2), newpage = F)
i = 4; print(SK_HFMD_sim, split = c(2,2,2,2), newpage = F)

dev.off()





#
#--- END ---
#




