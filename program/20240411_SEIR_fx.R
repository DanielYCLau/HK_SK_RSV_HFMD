




#
#--- PACKAGES ---
#





library("pomp")
library("coda")
library("lubridate")
library("scales")
library("doSNOW")
library("gridExtra")
library("reshape2")
library("pracma")
library("splines")
library("dplyr")

n.cores = parallel::detectCores() - 1

os = Sys.info()[["sysname"]]





#
#--- FUNCTIONS ---
#





pomp.dir = sprintf("~/%sDesktop/pomp_cfile/", ifelse(os == "Darwin", "", "../"))





pomp.v = as.numeric( substr(packageVersion("pomp"),1 ,3) )





win.fx = function(...) {
  if (os == "Windows") { windows(...) } else { quartz(...) }
}

file.fx = function(filename) {
  sprintf( "%s_%s", format(Sys.Date(), "%Y%m%d"), filename )
}





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



polygon.fx = function(x, LB, UB, ...) {
  # plot polygon
  polygon(x = c(x, rev(x)), y = c(LB, rev(UB)), border = NA, ...)
}





col.fx = function(ind, alpha = 1) {
  col = c(
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", 
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
  )
  return(alpha(col[ind], alpha))
}





SIR.scaling.fx = function(data, y, scale = T) {
  
  # calculate attack rate by season
  data = within(
    data, {
      year = year(WeekEnd)
      month = month(WeekEnd)
      season = ifelse(month < 10, year - 1, year)
    }
  )
  att.r = split(data, data[, "season"])
  
  att.r = sapply(att.r, function(x) { sum(x[, y]) })
  
  
  if (scale) { 
    # scaling 
    scale = min( 0.15 / att.r )
    # rate is around 15% - 50%  
    # att.r * scale
  } else {
    scale = 1
  }
  
  # output
  data = within(
    data, {
      scale       = scale
      scaled.data = get(y) * scale
    }
  )
  
  return(data)
}





RSV.age.scaling.fx = function(RSV.data, age.scale = NA) {
  
  #-- scaling for age --
  
  # only % RSV positive under 18 was reported since 22-Feb-2020 to 16-July-2022
  
  # According to Age-specific epidemic waves of RSV in Hong Kong from 2004 to 2013
  # (Yang. L, et al. Sci Rep.2015)
  # (https://www.nature.com/articles/srep10390/tables/1)
  # the % RSV positive under 18 was around 7.6% while for all ages was 4.2%
  
  # we thus scaled the % RSV positive under 18 to that for all ages for data since 2014
  
  age.0017.pos = (3322+177)/(33151+12814)
  age.all.pos  = (5012/120571)
  if (is.na(age.scale)) {
    age.scale    = age.all.pos / age.0017.pos
  }
  
  RSV.data = within(
    RSV.data, {
      p.RSV.all.age = p.RSV * ifelse(Age == "<18", age.scale, 1)
      RSV.proxy.GOPC = GOPC * p.RSV.all.age
      RSV.proxy.PMP  = PMP  * p.RSV.all.age
    }
  )  
  
  return(RSV.data)
}





data.formatting.fx = function(place, virus, scale) {
  
  
  
  # data loading
  if (place == "HK") {
    load(sprintf("../data/20240318_%s_%s_data.Rdata", place, virus))
  } else {
    load(sprintf("../data/20240402_%s_%s_data.Rdata", place, virus))
  }
  load(sprintf("../data/20240307_%s_demo_%s.Rdata", place, virus))
  load("../data/20240409_CHI.Rdata")
  index = index.list[[place]]
  index[, "index"] = zoo::na.locf(index[, "index"])
  assign("index", index, envir = .GlobalEnv)
  
  
  
  # data formatting
  
  if (virus == "HFMD") {
    
    if (place == "HK") {
      
      HFMD.data = SIR.scaling.fx(HFMD.data, "PMP", scale)
      
    } else if (place == "SK") {
      
      # HFMD.data = subset(HFMD.data, WeekEnd <= as.Date("2023-03-01"))
      HFMD.data = SIR.scaling.fx(HFMD.data, "HFMD.proxy", scale)
      
    }
    
    assign("HFMD.data", HFMD.data, envir = .GlobalEnv)
    
  } else if (virus == "RSV") {
    
    if (place == "HK") {
      
      RSV.data = RSV.age.scaling.fx(RSV.data)
      RSV.data = SIR.scaling.fx(RSV.data, "RSV.proxy.PMP", scale)
      
    } else if (place == "SK") {
      
      # RSV.data = subset(RSV.data, WeekEnd <= as.Date("2023-03-01"))
      RSV.data = SIR.scaling.fx(RSV.data, "RSV.proxy", scale)
      
    }
    
    assign("RSV.data", RSV.data, envir = .GlobalEnv)
    
  }
  
  assign("demo.sp.df", demo.sp.df, envir = .GlobalEnv)
  
  
}





global.C.fx = function(place, virus, nseas, nCOVbs) {
  
  # virus specific parameters
  if (virus == "RSV") {
    # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090094
    lat = 5/7 # rep(5/7, nAge)
    rec = 5/7 # rep(5/7, nAge)
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3481264/
    wan = 52.14 # rep( 52.14, nAge)
    
    # https://onlinelibrary.wiley.com/doi/full/10.1111/irv.12730
    # https://respiratory-research.biomedcentral.com/articles/10.1186/s12931-020-01456-3
    
    R0.UB = 5 # 10
    
  } else if (virus == "HFMD") {
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3987023/#:~:text=HFMD%20is%20a%20highly%20contagious,types%20of%20enterovirus%20%5B2%5D.
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418981/
    # https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-019-4153-6
    # https://www.cdc.gov/hand-foot-mouth/about/transmission.html
    lat = 5/7 # rep(   5/7, nAge)
    rec = 7/7 # rep(   7/7, nAge)
    wan = 52.14 # rep( 52.14, nAge) # assumed
    
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5117511/
    
    R0.UB = 10
    
  }
  
  
  
  rate.C = sprintf(
    "
    double lat = %f;
    double rec = %f;
    double wan = %f;
  ",
  lat, rec, wan
  )
  
  assign(    "lat",     lat, envir = .GlobalEnv)
  assign(    "rec",     rec, envir = .GlobalEnv)
  assign(    "wan",     wan, envir = .GlobalEnv)
  
  
  
  global.C = sprintf(
    "
    int    nseas  = %d; 
    int    nCOVbs = %d;
    double R0_UB  = %f;
    int    beta_ind = %d;
    int    b_c1_ind = %d;
    int    SKHFMD_ind = %d;
    ", 
    nseas, nCOVbs, R0.UB, 
    ifelse(place == "HK" & virus == "HFMD", 1, 0),
    ifelse(place == "HK" & virus == "RSV",  1, 0),
    ifelse(place == "SK" & virus == "HFMD",  1, 0)
  )
  
  
  

  
  # output
  # global.pomp = paste(global.C, age.C, CM.C, collapse = " ")
  global.pomp = paste(global.C, rate.C, collapse = " ")
  global.pomp = Csnippet(global.pomp)
  return(global.pomp)
}



# model for observed data
rMeas = Csnippet(
  "
  double Odis1 = exp( log_psi_1 );
  double Odis2 = exp( log_psi_2 );
  double Odis;  

  if (COVbs1 == 1.0) {
   Odis = Odis1;
  } else {
   Odis = Odis2;
  }

  obs  = rnbinom_mu(Odis, obs_mu);
"
)





# likelihood
dMeas = Csnippet(
  "
  double Odis1 = exp( log_psi_1 );
  double Odis2 = exp( log_psi_2 );
  double Odis;  

  if (COVbs1 == 1.0) {
   Odis = Odis1;
  } else {
   Odis = Odis2;
  }
  
  double obs_loglik;

  if ( ISNA( obs ) ) {
    obs_loglik = (give_log) ? 0.0 : 1.0 ;
  } else {
    obs_loglik = dnbinom_mu(obs, Odis, obs_mu, give_log);
  } 

  lik = obs_loglik ;

")






SEIR.model = Csnippet(
  "

  // Births
  double births = rpois(dt * all * birth);



  // rates
  double sigma = 1.0 / lat;    // rate for latent period
  double gamma = 1.0 / rec;    // recovery rate
  double rho   = 1.0 / wan;    // waning rate for infection-acquired immunity


  
  // transmissibility
  const double *beta_d = &beta1;
  double zeta[nseas];
  for (int j = 0; j < nseas; j++) {
    zeta[j] = beta_d[j];
  }
  zeta[0] = (beta_ind) ? (exp(zeta[0]) + 0.2 ) : zeta[0] ;
  zeta[1] = (SKHFMD_ind) ? (-exp(zeta[1])) : zeta[1];
  
  double r0     = R0_UB * expit(dot_product(nseas, &seas1, zeta));
  double beta_t = r0 * (gamma + death) * (sigma + death) / sigma;
  double seed;
  double covid_p;
  double lambda_t;

  
  
  // reduction in transmissibility
  double b_c_trans[nCOVbs];
  double b_c1_v2;
  b_c1_v2 = (b_c1_ind)? -exp(b_c1) : b_c1;

  
  b_c_trans[0] = 1.0;
  b_c_trans[1] = expit(b_c1_v2 - exp(b_c2) * index );
  b_c_trans[2] = expit(b_c1_v2 - exp(b_c2) * index + exp(b_c3));

  covid_p = dot_product(nCOVbs, &COVbs1, b_c_trans);

  if (COVbs1 == 1.0) {
    seed = 1.0e-4;
  } else {
    seed = 0.005 * expit(logit_seed);
  }

  lambda_t = beta_t * covid_p * fmax( I/pop, seed ); // important to have seed



  // transition between compartments

  // From class S
  double rateS[2], transS[2];
  rateS[0] = lambda_t;
  rateS[1] = death;
  
  reulermultinom(2, S, &rateS[0], dt, &transS[0]);
  
  S   += births;
  S   += - transS[0] - transS[1];
  E   +=   transS[0]; 
  Inc +=   transS[0];



  // From class E
  double rateE[2], transE[2];
  rateE[0] = sigma;
  rateE[1] = death;
    
  reulermultinom(2, E, &rateE[0], dt, &transE[0]);
  
  E += - transE[0] - transE[1];
  I +=   transE[0];



  // From class I
  double rateI[2], transI[2];
  rateI[0] = gamma;
  rateI[1] = death;
    
  reulermultinom(2, I, &rateI[0], dt, &transI[0]);

  I += - transI[0] - transI[1];
  R +=   transI[0];



  // From class R
  double rateR[2], transR[2];
  rateR[0] = rho;
  rateR[1] = death;
  
  reulermultinom(2, R, &rateR[0], dt, &transR[0]);
  
  R += - transR[0] - transR[1];
  S +=   transR[0]; 



  // From class Inc
  double scaling = expit( phi );
  obs_mu = Inc * scaling;



  // output state
  
  double tmp_N = S + E + I + R;

  S = round(pop * S / tmp_N);
  E = round(pop * E / tmp_N);
  I = round(pop * I / tmp_N);
  R = round(pop * R / tmp_N);

  ss   = S / pop; s_all = S / all;
  ee   = E / pop;
  ii   = I / pop;
  rr   = R / pop;
  iinc = Inc / pop;

  BBeta_t   = beta_t;
  LLambda_t = lambda_t;
  CCovid_p  = covid_p;

  RR0       = r0;
  RRe       = r0 * ss * covid_p;

"
)





rInit = Csnippet(
  "

  double Z[4];
  Z[0] = 0.7;
  Z[1] = 0.0; 
  Z[2] = 1.0e-4;
  Z[3] = 1.0 - Z[0] - Z[1] - Z[2];
  
  S   = nearbyint(pop * Z[0]);
  E   = nearbyint(pop * Z[1]);
  I   = nearbyint(pop * Z[2]);
  R   = nearbyint(pop * Z[3]);
  Inc = 0.0;
  
"
)




model.spec.fx = function(version, place, virus, stdt, eddt) { #version, place) {
  
  
  
  
  
  #-- data for model --
  
  if (virus == "HFMD") {
    data = HFMD.data
  } else if (virus == "RSV") {
    data = RSV.data
  }
  data = subset(data, WeekEnd >= as.Date(stdt) & WeekEnd <= as.Date(eddt))
  data = merge(x = data, y = demo.sp.df, by = "WeekEnd", all.x = T)
  
  pomp.data = with(
    data, {
      data.frame(
        "date" = WeekEnd,
        "time" = 1:nrow(data),
        "obs"  = round( scaled.data * pop ),
        "pop"  = pop,
        "all"  = all
      )
    }
  )
  pomp.data = transform(
    pomp.data,
    obs = ifelse(obs == 0, min(obs[obs > 0], na.rm = T), obs)
  )
  
  
  
  
  warm.up = 10 * 52
  warm.up.data = data.frame(
    "date" = rev( seq(data[1, "WeekEnd"] - 7, by = "-1 week", length.out = warm.up+1) ),
    "time" = -warm.up:0,
    "obs"  = rep(NA, warm.up+1),
    "pop"  = rep(NA, warm.up+1),
    "all"  = rep(NA, warm.up+1)
  )
  
  pomp.data = rbind(warm.up.data, pomp.data)
  assign("pomp.data", pomp.data, envir = .GlobalEnv)
  
  
  
  
  
  #-- time covariates --
  
  time.df.date = seq(pomp.data[1, "date"]-7, as.Date(eddt), by = "week")
  time.df      = data.frame( "time" = seq(-warm.up-1, by = 1, length.out = length(time.df.date) ) )
  
  
  
  # change in demography
  demo.df = subset(demo.sp.df, WeekEnd %in% time.df.date)
  
  if ( place == "HK" & virus == "RSV" ) {
    seas.df = data.frame(
      "seas1" = 1,
      "seas2" = sin( 2 * pi * time.df[, "time"] / 52.14  ),
      "seas3" = cos( 2 * pi * time.df[, "time"] / 52.14  ),
      "seas4" = sin( 2 * pi * time.df[, "time"] / 52.14 * 2 ),
      "seas5" = cos( 2 * pi * time.df[, "time"] / 52.14 * 2 )
    )
  } else { 
    seas.df = data.frame(
      "seas1" = 1,
      "seas2" = sin( 2 * pi * time.df[, "time"] / 52.14  ),
      "seas3" = cos( 2 * pi * time.df[, "time"] / 52.14  )
    )
  }
  assign("nSeas", ncol(seas.df), envir = .GlobalEnv)
  
  
  
  # COVID stringency index
  # if (place == "HK") {
  #   covid.t = data.frame(
  #     "date"    = time.df.date,
  #     "pre"     = 1 * ( time.df.date <= as.Date("2019-12-31") ),
  #     "covid"   = 1 * ( time.df.date  > as.Date("2019-12-31") & time.df.date  <= as.Date("2023-02-28")),
  #     "post"    = 1 * ( time.df.date  > as.Date("2023-02-28") )
  #   )
  # } else if (place == "SK") {
    
    
    post.w = seq(1, 0, length.out = sum( time.df.date  > (as.Date("2022-12-31") - 7) ) )
    
    covid.t = data.frame(
      "date"    = time.df.date,
      "pre"     = 1 * ( time.df.date <= as.Date("2019-12-31") ),
      "covid"   = 1 * ( time.df.date  > as.Date("2019-12-31") & time.df.date  <= as.Date("2022-12-31")),
      "post1"   = c(rep(0, sum( time.df.date  <= as.Date("2022-12-31") - 7 )), post.w),
      "post2"   = c(rep(0, sum( time.df.date  <= as.Date("2022-12-31") - 7 )), rev(post.w))
    )
    
    covid.t = with(
      covid.t,
      data.frame(
        "date"  = date,
        "pre"   = pre,
        "covid" = apply(cbind(covid, post1), 1, max),
        "post"  = post2
      )
    )
    
  # }

  
  
  
  
  

  covid.t = merge(
    x     = covid.t,
    y     = index,
    by.x  = "date",
    by.y  = "weekend",
    all.x = T
  )
  
  covid.t[, "index"] = with(covid.t, ifelse(is.na(index), 0, index) )

  covid.t = covid.t[, c("pre", "covid", "post", "index")]
  colnames(covid.t) = c("COVbs1", "COVbs2", "COVbs3", "index")

  assign("nCOVbs", sum(grepl("COVbs", colnames(covid.t))), envir = .GlobalEnv)
  
  
  
  
  # merge all covariate tgt
  
  cov.mat  = cbind(time.df, demo.df, seas.df, covid.t)
  cov.mat  = cov.mat[, !duplicated(colnames(cov.mat))]

  pomp.cov = covariate_table(cov.mat, order = "linear", times = "time")
  assign("cov.mat", cov.mat, envir = .GlobalEnv)
  
  
  
  
  
  #-- global variables --
  
  global.pomp = global.C.fx(
    place    = place,
    virus    = virus,
    nseas    = nSeas,
    nCOVbs   = nCOVbs
  )
  
  
  
  
  
  #-- model compartment and time series --
  
  state.names = c( # N is only allowed in stochastic
    "S", "E", "I", "R", "Inc",
    "ss", "ee", "ii", "rr", "iinc", "s_all",
    "BBeta_t", "LLambda_t", "CCovid_p", "obs_mu", 
    "RR0", "RRe"
  )
  
  obs.names   = c("obs")
  zero.names  = c("Inc")
  plot.names  = c(
    "ss", "ee", "ii", "rr", "iinc", "s_all",
    "obs_mu", "obs", 
    "BBeta_t", "LLambda_t", "CCovid_p", 
    "RR0", "RRe"
  )
  
  assign("state.names", state.names, envir = .GlobalEnv)
  assign("obs.names",   obs.names,   envir = .GlobalEnv)
  assign("plot.names",  plot.names,  envir = .GlobalEnv)
  
  
  
  
  
  #-- parameters --
  
  par.range = lapply(1:ncol(seas.df), function(x) { ifelse(x == 1, 2, 1) * c(-1, 1) } )
  names(par.range) = sprintf("beta%d", 1:ncol(seas.df))
  par.range = c(
    par.range, 
    list(
      "log_psi_1"  = c(-2, 2),
      "phi"        = c(-3,-1),
      "b_c1"       = c(-2, 2),
      "b_c2"       = c(-10,-5),
      "logit_seed" = c(-5,-1),
      "log_psi_2"  = c(-2, 2),
      "b_c3"       = c(-2, 2)
    )
  )
  # if (place == "HK") {
  #   par.range = c( par.range, list("b_c3" = c(-2, 2)) )
  # }

  assign("par.range", par.range, envir = .GlobalEnv)
  
  par.names = names(par.range)
  assign("par.names", par.names, envir = .GlobalEnv)
  
  
  
  
  
  # par.rw.sd for mif2
  t.2020.ind = subset(cov.mat,  COVbs1  < 1)[1, "time"]
  # if (place == "HK") { t.2023.ind = subset(cov.mat,  COVbs3 == 1)[1, "time"] }
  t.2023.ind = subset(cov.mat,  COVbs3 > 0)[1, "time"] 
  
  par.rw.sd.t = sapply(
    par.names,
    function(x) {
      if (x %in% c("b_c1", "b_c2", "log_psi_2", "logit_seed")) {
        t = t.2020.ind
      } else if (x == "b_c3") {
        t = t.2023.ind
      } else {
        t = 0
      }
      return(t)
    }
  )
  par.rw.sd.text = sprintf(
    "rw%ssd(%s)",
    ifelse(pomp.v < 5, ".", "_"),
    paste( sprintf( "%s = ifelse(time >= %s, 0.02, 0)", par.names, par.rw.sd.t), collapse = ", " )
  )
  par.rw.sd      = eval(parse(text = par.rw.sd.text))
  assign("par.rw.sd", par.rw.sd, envir = .GlobalEnv)
  
  
  
  
  
  # Prior for mcmc
  dPri = Csnippet(
    sprintf(
      "
      lik = %s;
      lik = (give_log) ? lik : exp(lik);
      ",
      paste( sprintf("dunif(%s, -10, 10, 1)", par.names), collapse = " + " )
    )
  )
  # }
  assign("dPri", dPri, envir = .GlobalEnv)
  
  
  
  
  
  
  
  #-- model --
  
  pomp.euler = euler(step.fun = SEIR.model, delta.t = 1/7)
  
  pomp.model = pomp(
    cdir       = pomp.dir,
    cfile      = sprintf("%s_%s_SEIR_%s", version, place, virus),  
    data       = pomp.data[, c("time", "obs")],
    times      = "time",
    t0         = cov.mat[1, "time"],
    globals    = global.pomp,
    rprocess   = pomp.euler,
    rmeasure   = rMeas,
    dmeasure   = dMeas,
    dprior     = dPri,
    covar      = pomp.cov,
    rinit      = rInit,
    statenames = state.names,
    obsnames   = obs.names, 
    accumvars  = zero.names, 
    paramnames = par.names
  )
  
  return(pomp.model)
}





proposal.fx = function(...) {
  if (pomp.v < 5) {
    return(mvn.rw.adaptive(...))
  } else {
    return(mvn_rw_adaptive(...))    
  }
}





plot.fx = function(state, sim, pre = F) {
  
  # state = "R0"; sim = if.sim; pre = F
  time = sim[, "time"]
  
  if ( any(colnames(sim) %in% state) ) {
    col.names = state
  } else {
    col.ind   = grep(sprintf("^%s_", state), colnames(sim))
    col.names = colnames(sim)[col.ind]
  } 
  
  
  if ( length(col.names) > 1 ) {
    sim[, col.names] = sim[, col.names]
    sim.sum = aggregate(sim[, col.names], by = list("time" = time), FUN = mean)
    sim.sum = within(sim.sum, { if (!state %in% c("R0", "Re")) { All = rowSums(sim.sum[, col.names]) } } )
  } else {
    sim.sum = aggregate(sim[, col.names], by = list("time" = time), FUN = CI95.fx)
    sim.sum = do.call(data.frame, sim.sum)
  }
  
  if (!pre) { sim.sum = subset(sim.sum, time > 0)}
  
  
  time  = merge(x = sim.sum, y = pomp.data, by = "time", all.x = T)
  time  = decimal_date(time[, "date"])
  sim.x = subset(sim.sum, select = -time)
  
  
  with(
    sim.x, {
      
      if ( state == "obs" ) {
        plot(
          NULL, 
          xlim = range(time), 
          ylim = c(0, max(pomp.data[, "obs"] * 2, na.rm = T)), 
          main = state,
          ylab = "",
          xaxt = "n"
        )
      } else {
        plot(
          NULL, 
          xlim = range(time), 
          ylim = c(0, max(sim.x, na.rm = T)), 
          main = state,
          ylab = "",
          xaxt = "n"
        )
      }
      
      axis.t = round(min(time)):2050
      abline(v = axis.t, col = alpha(1, 0.1))
      axis(1, at = axis.t, labels = F)
      axis(1, at = axis.t + 0.5, labels = axis.t, tick = F)
      
      if (length(col.names) == 1) { 
        lines(x = time, y = x.mean)
        polygon.fx(x = time, LB = x.LB, UB = x.UB, col = alpha(1, 0.1))
      } else {
        for (i in 1:ncol(sim.x)) {
          if ( grepl("All|all", colnames(sim.x)[i]) ) {
            lines(x = time, y = sim.x[,i], col = 1, lwd = 2)
          } else {
            lines(x = time, y = sim.x[,i], col = i)
          }
        }  
      }
    }
  )
  if ( state %in% c("R0", "Re") ) { abline(h = 1, col = 2) }
  if (pre) { abline(v = 0, col = 2) }
  
}




traceplot.fx = function(fit, par, alpha) {
  
  # fit = if.fit
  traces = lapply(fit, function(x) { x[["fit"]]@traces[, par] } )
  traces = do.call(rbind, traces)
  
  Npar  = nrow(traces)
  iter  = ncol(traces)
  y.lim = range(traces)
  
  plot(NULL, xlim = c(1, iter), ylim = y.lim, main = par, ylab = "")
  
  a = c(0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999)
  sapply(
    a, function(a.i) { 
      LB.i = apply(traces, 2, quantile, a.i / 2)
      UB.i = apply(traces, 2, quantile, 1 - a.i / 2)
      
      polygon.fx(x = 1:iter, LB = LB.i, UB = UB.i, col = alpha(1, alpha))
    }
  )
}





# profile plot for IF2
profile.plot.fx = function(fit, par, thres) {
  
  est = lapply(fit, function(x) { x[["est"]] } )
  est = do.call(rbind, est)
  
  max.loglik = max(est[, "loglik"], na.rm = T)
  if.est     = subset(est, loglik == max.loglik)
  est        = subset(est, loglik >= max.loglik - thres)
  
  profile.thres = qchisq(0.95, df = 1)/2
  
  x.lim = range(est[, par])
  
  with(est, plot(x = get(par), y = loglik, main = par, ylab = "") )
  with(
    if.est, {
      points(x = get(par), y = loglik, col = 2, pch = 16)
      abline(v = get(par), h = loglik, col = 2, lty = 2)
      abline(h = max.loglik - profile.thres, col = 3, lty = 2)
    }
  )
  
}




pairs.fx = function(fit, thres = 5) {
  
  est = lapply(fit, function(x) { x[["est"]] } )
  est = do.call(rbind, est)
  
  est = na.omit(est)
  est = subset(est, loglik >= max(loglik) - thres)
  pairs(est)
  
}





sim.fx = function(fit.list) {
  
  ind = which.max(sapply(fit.list, function(x) { x[["loglik"]] } ))
  fit.i = fit.list[[ind]]
  
  if2.model = pomp(
    data   = pomp.model,
    params = fit.i[["est"]]
  )
  assign("if2.model", if2.model, envir = .GlobalEnv)
  
  
  
  sim = simulate(
    object = if2.model,
    nsim   = 1000,
    seed   = 555,
    format = "data.frame"
  )
  
  return(sim)
}





plot.mcmc.fx = function(trace) {
  
  # trace = pmcmc.trace
  for (i in c("loglik", par.names) ) {
    traceplot(trace[, i], main = i)
    
    densplot(trace[, i], main = i, lwd = 2)
    
    abline(v = summary(trace)[["statistics"]][i, "Mean"], col = 2, lwd = 2)
    abline(v = summary(trace)[["quantiles"]][i, c("2.5%", "50%", "97.5%")], col = 2, lwd = 2, lty = 3)
    
    for (j in 1:length(trace)) {
      tmp.j = trace[[j]]
      lines(density(tmp.j[, i]), col = alpha(j, 0.25))
    }
    
  }
}




Neff.condlogLik.fx = function(fit, pre) {
  
  # fit = if.fit
  ind = which.max(sapply(fit, function(x) { x[["loglik"]][1] } ))
  fit = fit[[ind]][["fit"]]
  
  Neff        = eff_sample_size(fit)
  cond.loglik = cond_logLik(fit)
  
  tmp = within(
    pomp.data, {
      dec.date    = decimal_date(date)
      cond.loglik = cond.loglik
      Neff        = Neff
    }
  )
  
  if (!pre) { tmp = subset(tmp, time > 0)}
  
  plot(cond.loglik ~ dec.date, data = tmp, type = "l")
  abline(v = 2000:2030, col = alpha(1, 0.1))
  plot(Neff        ~ dec.date, data = tmp, type = "l")
  abline(v = 2000:2030, col = alpha(1, 0.1))
  
}





pmcmc.init.fx = function(fit.list) {
  
  # fit.list = if.fit
  pmcmc.init = lapply(fit.list, function(x) { x[["est"]] } )
  pmcmc.init = do.call(rbind, pmcmc.init)
  
  pmcmc.init = pmcmc.init[with(pmcmc.init, order(-loglik)), ]
  
  # if (place == "SK" & virus == "RSV") {
    pmcmc.init = subset(pmcmc.init, loglik.se < 1)
    init.ind   = with(pmcmc.init, loglik >= (max(loglik) - 100))
  # } else {
  #   pmcmc.init = subset(pmcmc.init, loglik.se < 0.1)
  #   init.ind   = with(pmcmc.init, loglik >= (max(loglik) - 5))
  # }
  
  if ( sum(init.ind) <= 10 ) {
    pmcmc.init = pmcmc.init[1:10, ]
  } else {
    pmcmc.init = pmcmc.init[init.ind, ]
  }
  pmcmc.init = subset(pmcmc.init, select = -c(loglik, loglik.se))
  for (k in colnames(pmcmc.init)) { # for dPri
    pmcmc.init[, k] = ifelse(pmcmc.init[, k] >  10,  9.9, pmcmc.init[, k])
    pmcmc.init[, k] = ifelse(pmcmc.init[, k] < -10, -9.9, pmcmc.init[, k])
  }

  return(pmcmc.init)
}





pmcmc.sum.fx = function(trace) {
  
  # trace = pmcmc.trace
  n    = sum( sapply(trace, nrow) )
  thin = attr(trace[[1]], "mcpar")[3]
  
  R    = gelman.diag(trace[, par.names], autoburnin = F, transform = T)[["psrf"]]
  R    = as.data.frame(R)
  R    = apply(R, 2, function(x) { sprintf("%.2f", x) } )
  colnames(R) = c("R", "R.UB")
  
  stat = summary(trace)[["statistics"]][par.names, c("Mean", "SD")]
  stat = as.data.frame(stat)
  stat = apply(stat, 2, function(x) { sprintf("%.2f", x) } )
  
  q    = summary(trace)[["quantiles"]][par.names, c("2.5%", "50%", "97.5%")]
  q    = as.data.frame(q)
  colnames(q) = c("LB", "Q50", "UB")
  q[, "CI95"] = with(q, sprintf("(%.2f, %.2f)", LB, UB)) 
  
  sum = cbind(
    "par"  = par.names, 
    "N"    = n, 
    "thin" = thin, 
    R, 
    stat, 
    "Q50"  = sprintf("%.2f", q[, "Q50"]), 
    "CI95" = q[, "CI95"]
  )
  sum = as.data.frame(sum)
  rownames(sum) = NULL
  
  return(sum)
}




pmcmc.sim.fx = function(trace) {
  
  # trace = pmcmc.trace
  # trace = pmcmc.coda
  if (is.data.frame(trace)) {
    coda = t(trace)
  }else if (is.list(trace)) {
    coda = do.call(rbind, trace)
    coda = t(coda)
  } 
  
  set.seed(555)
  n.sim = ifelse( ncol(coda) < 10000, ncol(coda), 10000)
  coda  = coda[, sample(1:ncol(coda), n.sim, replace = F)]
  colnames(coda) = 1:n.sim
  
  sim = simulate(
    object = pomp.model,
    params = coda,
    nsim   = 1,
    seed   = 555,
    format = "data.frame"
  )
  
  return(sim)
}





if (os == "Linux") { print( sessionInfo() ) }





#
#--- END ---
#




