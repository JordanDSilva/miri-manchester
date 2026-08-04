library(celestial)
library(ProSpect)
library(data.table)
library(foreach)
library(stringr)
library(magicaxis)
library(Highlander)
library(ks)
library(MASS)
library(fields)
library(lamW)
library(rhdf5)
library(hyper.fit)
library(rgl)
library(matrixStats)
library(LaplacesDemon)
library(mvtnorm)  # for dmvnorm

prospect_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

SMBH_df = data.frame(fread("~/Documents/Epochs/data/prospect/prospect_robust_SMBH_results.csv"))
BHARSFR_df = data.frame(fread("~/Documents/Epochs/data/prospect/prospect_robust_BHAR_SFR_ratio.csv"))

LAGN2BHAR = ( 3.154e7 * 1e-7 * (0.1 * 3e8^2)**-1 / 2e30)

lbt_bins = seq(8.0, 13.0, 0.5)
zbins = approx(
  x = cosdistTravelTime(z = seq(0, 40, 0.01), ref = "737"),
  y = seq(0, 40, 0.01),
  xout = lbt_bins
)$y
zmids = zbins[-1] - diff(zbins)*0.5

sample_idx = prospect_cat['zQ50'] >= min(zbins) & prospect_cat['zQ50'] < max(zbins) & !is.na(prospect_cat["zQ50"])
quenched_idx = prospect_cat["SFR10Q50"][sample_idx]/prospect_cat["StellarMassQ50"][sample_idx] <= 0.2/(1e9*cosdistUniAgeAtz(prospect_cat["zQ50"][sample_idx], ref = "737")) # probably just use this for upper limits on the SFR10 plots 

prospect_sample = prospect_cat[sample_idx, ]
BHARSFR_df_sample = BHARSFR_df[sample_idx, ]

dbl_pl = function(x, parm){
  M0    = parm[1]
  alpha = parm[2]
  beta  = parm[3]
  S0    = parm[4]
  
  xx = 10^(x - M0)
  
  y = S0 - log10(
    (xx)^(-1*alpha) + (xx)^(-1*beta)
  )
  return(y)
}

PGF = function(Data){
  M0    = rnorm(1, 7, 0.5)
  alpha = rnorm(1, 1, 0.3)
  log_delta = rnorm(1, log(0.5), 0.5)
  # delta = rgamma(1, shape = 2, rate = 2)   # ensures positivity
  S0    = rnorm(1, 1, 0.3)
  xtrue = rnorm(Data$N, Data$x_obs, Data$xerr)
  return(c(M0, alpha, log_delta, S0, xtrue))
}

Model = function(parm, Data){
  
  M0    = parm[1]
  alpha = parm[2]
  log_delta = parm[3]
  S0    = parm[4]
  xtrue = parm[5:(4 + Data$N)]
  
  delta = exp(log_delta)
  beta = alpha - delta
  
  # if (alpha <= beta){
  #   return(list(LP = -Inf, Dev = Inf, Monitor = c(-Inf, M0, alpha, beta, S0),
  #               yhat = rep(NA, Data$N), parm = parm))
  # }
  
  # --- Priors on model parameters ---
  M0.prior    = dnorm(M0,    7, 2, log = TRUE)
  alpha.prior = dnorm(alpha, 1, 2, log = TRUE)
  log_delta.prior = dnorm(log_delta, log(0.5), 1, log = TRUE)
  # delta.prior = dgamma(delta, shape = 2, rate = 2, log = TRUE)  # keeps delta > 0
  # beta.prior  = dnorm(beta,  1, 2, log = TRUE)
  S0.prior    = dnorm(S0,    1, 2, log = TRUE)
  
  # --- Model prediction at xtrue ---
  y_model = dbl_pl(xtrue, c(M0, alpha, beta, S0))
  
  # --- 2D Gaussian likelihood per point ---
  # LL = 0
  # for (j in seq_len(Data$N)){
  #   obs_vec   = c(Data$x_obs[j], Data$y_obs[j])
  #   model_vec = c(xtrue[j],      y_model[j])
  #   Sigma_j   = Data$covarray[, , j]
  #   
  #   LL = LL + dmvnorm(obs_vec, mean = model_vec, sigma = Sigma_j, log = TRUE)
  # }
  # --- Extract covariance components as vectors (length N) ---
  sx2 = Data$covarray[1, 1, ]   # xerr^2 per point
  sy2 = Data$covarray[2, 2, ]   # yerr^2 per point
  sxy = Data$covarray[1, 2, ]   # covariance per point (0 if uncorrelated)
  
  dx = Data$x_obs - xtrue
  dy = Data$y_obs - y_model
  
  det_Sigma = sx2 * sy2 - sxy^2
  
  # --- Closed-form bivariate normal log-density, vectorized ---
  quadform = (dx^2 * sy2 - 2 * dx * dy * sxy + dy^2 * sx2) / det_Sigma
  
  LL_vec = -log(2 * pi) - 0.5 * log(det_Sigma) - 0.5 * quadform
  
  LL = sum(LL_vec)
  
  LP = LL + M0.prior + alpha.prior + log_delta.prior + S0.prior
  
  Modelout = list(
    LP = LP,
    Dev = -2 * LL,
    Monitor = c(LP, M0, alpha, beta, delta, S0),
    yhat = y_model,
    parm = parm
  )
  return(Modelout)
}

# --- Root-finding function for a single parameter draw ---
find_zero = function(parm, interval = c(1, 15)){
  f = function(x) dbl_pl(x, parm)
  
  # check that a root actually exists in this interval (opposite signs at endpoints)
  f_lo = f(interval[1])
  f_hi = f(interval[2])
  
  if (is.na(f_lo) || is.na(f_hi) || sign(f_lo) == sign(f_hi)){
    return(NA)
  }
  
  root = uniroot(f, interval = interval)$root
  return(root)
}
find_zero_from_grid = function(x_grid, y_vals){
  # locate sign change
  sign_change = which(diff(sign(y_vals)) != 0)
  
  if (length(sign_change) == 0){
    return(NA)
  }
  
  # use the first crossing (adjust if multiple crossings are possible/expected)
  idx = sign_change[1]
  
  # linear interpolation between the bracketing grid points
  x_zero = approx(
    x = y_vals[idx:(idx+1)],
    y = x_grid[idx:(idx+1)],
    xout = 0
  )$y
  
  return(x_zero)
}

x_grid = seq(-15, 15, 0.01)

df_fit = list()
for(i in 1:length(zmids)){
  
  zidx = prospect_sample$zQ50 >= zbins[i] & prospect_sample$zQ50 < zbins[i + 1] & BHARSFR_df_sample$lowerlim == FALSE & prospect_sample$FitType == "StellarAGN"

  # --- Data ---
  x_obs = log10(prospect_sample$StellarMassQ50[zidx])
  y_obs = log10(BHARSFR_df_sample$BHARSFR10Q50[zidx])
  
  xerr = (log10(prospect_sample$StellarMassQ84[zidx]) - log10(prospect_sample$StellarMassQ16[zidx])) * 0.5
  yerr = (log10(BHARSFR_df_sample$BHARSFR10Q84[zidx]) - log10(BHARSFR_df_sample$BHARSFR10Q16[zidx])) * 0.5
  
  N = length(x_obs)
  
  covarray = array(0, dim = c(2, 2, N))
  for (j in 1:N){
    read_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/",
      "/galID_", prospect_sample$CATID[zidx][j], ".rds"
    )
    read_data = readRDS(read_stub)
    post = read_data$post[, c("StellarMass", "SFR10", "AGNlum")]
    
    x = as.numeric(log10(post[, "StellarMass"]))
    y = as.numeric(log10(post[, "SFR10"] / (LAGN2BHAR * post[, "AGNlum"])))
    
    ok = is.finite(x) & is.finite(y)
    
    if (sum(ok) < 2){
      # not enough valid samples to compute covariance — fall back to something safe
      warning(sprintf("j = %d: insufficient finite samples (%d), using diagonal fallback", j, sum(ok)))
      covarray[, , j] = diag(c(var(x[is.finite(x)], na.rm = TRUE),
                               var(y[is.finite(y)], na.rm = TRUE)))
      next
    }
    
    covarray[, , j] = cov(cbind(x[ok], y[ok]))
  }
  
  # 
  # covarray = array(0, dim = c(2, 2, N))
  # for (j in seq_len(N)){
  #   covarray[, , j] = matrix(
  #     c(xerr[j]^2, 0,
  #       0, yerr[j]^2),
  #     nrow = 2
  #   )
  # }
  
  Data = list(
    x_obs = x_obs,
    y_obs = y_obs,
    xerr  = xerr,
    yerr  = yerr,
    N     = N,
    covarray = covarray,
    # mon.names  = c("LP", "M0", "alpha", "beta", "S0"),
    # parm.names = c("M0", "alpha", "beta", "S0", paste0("xtrue", 1:N)),
    PGF   = PGF
  )
  Data$mon.names  = c("LP", "M0", "alpha", "beta", "delta", "S0")
  Data$parm.names = c("M0", "alpha", "log_delta", "S0", paste0("xtrue", 1:N))
  
  Initial.Values = PGF(Data)
  
  Fit = LaplacesDemon(
    Model = Model,
    Data  = Data,
    Initial.Values = Initial.Values,
    Covar = NULL,
    Iterations = 20000,
    Status = 1000,
    Thinning = 10,
    Algorithm = "CHARM",       # adaptive algorithm, good default for this kind of problem
    Specs = NULL
  )
  
  # --- Extract posterior on model parameters only (drop latent x_true columns) ---
  post = Fit$Monitor[, c("M0", "alpha", "beta", "S0")]
  
  # --- Grid over which to evaluate the fit ---
  y_samples = apply(post, 1, function(parm) dbl_pl(x_grid, parm))
  # y_samples is a matrix: nrow = length(x_grid), ncol = n_posterior_samples
  
  # --- Compute quantiles at each x point across posterior samples ---
  y_quantiles = t(apply(y_samples, 1, quantile, probs = c(0.50, 0.16, 0.84), na.rm = TRUE))
  colnames(y_quantiles) = c("Q50", "Q16", "Q84")
  
  df_quantiles_post = colQuantiles(
    x = post,
    probs = c(0.5, 0.16, 0.84),
    na.rm = TRUE
  )
  
  # --- Apply across all posterior samples ---
  Mbreak_quantiles = apply(y_quantiles, 2, function(x){find_zero_from_grid(x_grid, y_vals = x)})
  Mbreak_quantiles[c("Q16", "Q84")] = Mbreak_quantiles[c("Q84", "Q16")]
  Mbreak_quantiles[is.na(Mbreak_quantiles)] = -999
  
  df_quantiles = rbind(df_quantiles_post, "Mbreak" = Mbreak_quantiles)
  colnames(df_quantiles) = c("Q50", "Q16", "Q84")
  
  df_lo = df_quantiles[, "Q50"] - df_quantiles[, "Q16"]
  df_hi = df_quantiles[, "Q84"] - df_quantiles[, "Q50"]
  
  df = c(
    setNames(df_quantiles[, "Q50"], paste0(rownames(df_quantiles), "Q50")),
    setNames(df_quantiles[, "Q16"], paste0(rownames(df_quantiles), "Q16")),
    setNames(df_quantiles[, "Q84"], paste0(rownames(df_quantiles), "Q84")),
    setNames(df_lo,                 paste0(rownames(df_quantiles), "_lo")),
    setNames(df_hi,                 paste0(rownames(df_quantiles), "_hi"))
  )

  y_quantiles = cbind(y_quantiles, "x" = x_grid)
  
  # --- Base plot ---
  pdf(
    paste0(
      "~/Documents/Epochs/plots/fit_bhar_sfr/z", i, ".pdf"
    ),
    width = 8,
    height = 6
  )
  magplot(
    log10(prospect_sample$StellarMassQ50[zidx]),
    log10(BHARSFR_df_sample$BHARSFR10Q50[zidx]),
    xlim = c(1, 12),
    ylim = c(-5, 5),
    xlab = "log10(Stellar Mass)",
    ylab = "log10(BHAR/SFR)"
  )
  magerr(
    log10(prospect_sample$StellarMassQ50[zidx]),
    log10(BHARSFR_df_sample$BHARSFR10Q50[zidx]),
    ylo = (log10(BHARSFR_df_sample$BHARSFR10Q84[zidx]) - log10(BHARSFR_df_sample$BHARSFR10Q16[zidx])) * 0.5,
    xlo = (log10(prospect_sample$StellarMassQ84[zidx]) - log10(prospect_sample$StellarMassQ16[zidx])) * 0.5
  )
  legend(
    x = "topleft",
    legend = paste0(
      round(zbins[i],3), " < z <", round(zbins[i+1], 3)
    )
  )
  
  # --- Shaded 16th-84th percentile polygon for the fit ---
  polygon(
    x = c(x_grid, rev(x_grid)),
    y = c(y_quantiles[, "Q16"], rev(y_quantiles[, "Q84"])),
    col = rgb(0, 0, 1, alpha = 0.25),   # transparent blue
    border = NA
  )
  
  # --- Median fit line ---
  lines(
    x_grid,
    y_quantiles[, "Q50"],
    col = "blue",
    lwd = 2
  )
  
  # --- Mbreak ---
  points(
    df["MbreakQ50"],
    0,
    pch = 17,
    col = "red",
    cex = 1.5
  )
  magerr(
    df["MbreakQ50"],
    0,
    xlo = df["Mbreak_lo"],
    xhi = df["Mbreak_hi"],
    col = "red"
  )
  dev.off()
  
  df_fit[[i]] = list(
    "df" = df,
    "fit" = y_quantiles
  )
}

df_cat = data.frame(do.call(
  rbind, lapply(df_fit, function(x)x$df)
))
saveRDS(
  df_fit,
  "~/Documents/Epochs/data/prospect/prospect_robust_SMBH_results.rds"
)

fwrite(
  df_cat, 
  "~/Documents/Epochs/data/prospect/BHAR_SFR_dpl_fit_results.csv"
)

df_fit = readRDS("~/Documents/Epochs/data/prospect/prospect_robust_SMBH_results.rds")
df_cat = data.frame(fread("~/Documents/Epochs/data/prospect/BHAR_SFR_dpl_fit_results.csv"))

h5file <- path.expand("~/Documents/Epochs/data/analysis_data.h5")
# Create the file if it doesn't exist yet; otherwise we append into it
if (!file.exists(h5file)) {
  h5createFile(h5file)
}

h5delete(h5file, "SFRBHAR")
h5createGroup(h5file, "SFRBHAR")
h5write(df_cat, file = h5file, name = paste0("SFRBHAR/Params"))
h5createGroup(h5file, "SFRBHAR/fits")
for(i in 1:length(df_fit)){
  h5write(
    data.frame(df_fit[[i]]$fit), 
    file = h5file, 
    name = paste0("SFRBHAR/fits/", i)
  )
}

foo = h5read(
  h5file, "SFRBHAR/fits/1"
)

magplot(
  zmids,
  df_cat$MbreakQ50,
  ylim = c(-0.5, 8),
  xlim = c(0, 10),
  ylab = "Mbreak"
)
magerr(
  zmids,
  df_cat$MbreakQ50,
  ylo = df_cat$Mbreak_lo,
  yhi = df_cat$Mbreak_hi
)
abline(h = 7)

magplot(
  zmids,
  df_cat$M0Q50,
  ylim = c(-0.5, 12),
  xlim = c(0, 10),
  ylab = "M0"
)
magerr(
  zmids,
  df_cat$M0Q50,
  ylo = df_cat$M0_lo,
  yhi = df_cat$M0_hi
)

magplot(
  zmids,
  df_cat$alphaQ50,
  ylim = c(-2.5, 2.5),
  xlim = c(0, 10),
  ylab = "alpha"
)
magerr(
  zmids,
  df_cat$alphaQ50,
  ylo = df_cat$alpha_lo,
  yhi = df_cat$alpha_hi
)

magplot(
  zmids,
  df_cat$betaQ50,
  ylim = c(-2.5, 2),
  xlim = c(0, 10),
  ylab = "beta"
)
magerr(
  zmids,
  df_cat$betaQ50,
  ylo = df_cat$beta_lo,
  yhi = df_cat$beta_hi
)




# --- Parameter vector layout ---
# parm[1:4]   = M0, alpha, beta, S0        (model parameters)
# parm[5:(4+N)] = x_true_i                  (latent true x for each point)

# PGF = function(Data){
#   # Parameter Generating Function - initial values for MCMC
#   M0    = rnorm(1, 7, 0.5)
#   alpha = rnorm(1, 1, 0.3)
#   beta  = rnorm(1, 1, 0.3)
#   S0    = rnorm(1, 1, 0.3)
#   xtrue = rnorm(Data$N, Data$x_obs, Data$xerr)
#   return(c(M0, alpha, beta, S0, xtrue))
# }
# 
# Model = function(parm, Data){
#   
#   M0    = parm[1]
#   alpha = parm[2]
#   beta  = parm[3]
#   S0    = parm[4]
#   xtrue = parm[5:(4 + Data$N)]
#   
#   # --- Priors on model parameters (adjust ranges as appropriate) ---
#   M0.prior    = dnorm(M0,    7, 2,   log = TRUE)
#   alpha.prior = dnorm(alpha, 1, 1,   log = TRUE)
#   beta.prior  = dnorm(beta,  1, 1,   log = TRUE)
#   S0.prior    = dnorm(S0,    1, 2,   log = TRUE)
#   
#   # --- Prior on alpha/beta
#   beta.prior = 
#   
#   # --- Prior on latent x_true: encodes that we expect it near the observed x ---
#   # (this is effectively the "measurement model" linking x_true to x_obs via xerr)
#   xtrue.prior = sum(dnorm(xtrue, Data$x_obs, Data$xerr, log = TRUE))
#   
#   LP.prior = M0.prior + alpha.prior + beta.prior + S0.prior + xtrue.prior
#   
#   # --- Likelihood: y_obs given model evaluated at x_true ---
#   mu = dbl_pl(xtrue, c(M0, alpha, beta, S0))
#   LL = sum(dnorm(Data$y_obs, mean = mu, sd = Data$yerr, log = TRUE))
#   
#   LP = LL + LP.prior
#   
#   Modelout = list(
#     LP = LP,
#     Dev = -2 * LL,
#     Monitor = c(LP, M0, alpha, beta, S0),
#     yhat = mu,
#     parm = parm
#   )
#   return(Modelout)
# }

# --- Parameter vector layout now: M0, alpha, delta, S0, xtrue... ---
# where delta > 0 and beta = alpha - delta



# Model = function(parm, Data){
#   
#   M0    = parm[1]
#   alpha = parm[2]
#   delta = parm[3]
#   S0    = parm[4]
#   xtrue = parm[5:(4 + Data$N)]
#   
#   beta = alpha - delta
#   
#   # --- Priors ---
#   M0.prior    = dnorm(M0,    7, 2, log = TRUE)
#   alpha.prior = dnorm(alpha, 1, 2, log = TRUE)
#   delta.prior = dgamma(delta, shape = 2, rate = 2, log = TRUE)  # keeps delta > 0
#   S0.prior    = dnorm(S0,    1, 2, log = TRUE)
#   
#   xtrue.prior = sum(dnorm(xtrue, Data$x_obs, Data$xerr, log = TRUE))
#   
#   LP.prior = M0.prior + alpha.prior + delta.prior + S0.prior + xtrue.prior
#   
#   mu = dbl_pl(xtrue, c(M0, alpha, beta, S0))
#   LL = sum(dnorm(Data$y_obs, mean = mu, sd = Data$yerr, log = TRUE))
#   
#   LP = LL + LP.prior
#   
#   # if(is.infinite(LP)){
#   #   LP = -99999
#   # }
#   # 
#   Modelout = list(
#     LP = LP,
#     Dev = -2 * LL,
#     Monitor = c(LP, M0, alpha, beta, delta, S0),
#     yhat = mu,
#     parm = parm
#   )
#   return(Modelout)
# }