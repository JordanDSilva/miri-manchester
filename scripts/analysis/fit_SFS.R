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
prospect_sample = prospect_sample[!quenched_idx, ]

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
  M0    = rnorm(1, 10, 0.5)
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
  
  # --- Priors on model parameters ---
  M0.prior    = dnorm(M0,    10, 2, log = TRUE)
  alpha.prior = dnorm(alpha, 1, 2, log = TRUE)
  log_delta.prior = dnorm(log_delta, log(0.5), 1, log = TRUE)
  S0.prior    = dnorm(S0,    1, 2, log = TRUE)
  
  # --- Model prediction at xtrue ---
  y_model = dbl_pl(xtrue, c(M0, alpha, beta, S0))
  
  # --- 2D Gaussian likelihood per point ---
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


x_grid = seq(-15, 15, 0.01)

df_fit = list()
for(i in 1:length(zmids)){
  
  zidx = prospect_sample$zQ50 >= zbins[i] & prospect_sample$zQ50 < zbins[i + 1] 
  
  # --- Data ---
  x_obs = log10(prospect_sample$StellarMassQ50[zidx])
  y_obs = log10(prospect_sample$SFR10Q50[zidx])
  
  xerr = (log10(prospect_sample$StellarMassQ84[zidx]) - log10(prospect_sample$StellarMassQ16[zidx])) * 0.5
  yerr = (log10(prospect_sample$SFR10Q84[zidx]) - log10(prospect_sample$SFR10Q16[zidx])) * 0.5
  
  N = length(xerr)
  
  covarray = array(0, dim = c(2, 2, N))
  for (j in 1:N){
    read_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/",
      "/galID_", prospect_sample$CATID[zidx][j], ".rds"
    )
    read_data = readRDS(read_stub)
    post = read_data$post[, c("StellarMass", "SFR10")]
    
    x = as.numeric(log10(post[, "StellarMass"]))
    y = as.numeric(log10(post[, "SFR10"]))
    
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
  df_quantiles = df_quantiles_post
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
      "~/Documents/Epochs/plots/rough-plots/fit_sfs/z", i, ".pdf"
    ),
    width = 8,
    height = 6
  )
  magplot(
    x_obs,
    y_obs,
    xlim = c(1, 12),
    ylim = c(-5, 5),
    xlab = "log10(Stellar Mass)",
    ylab = "log10(SFR)"
  )
  magerr(
    x_obs,
    y_obs,
    xlo = xerr,
    ylo = yerr
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
  dev.off()
  
  df_fit[[i]] = list(
    "df" = df,
    "fit" = y_quantiles
  )
}
 
# foo = h5read(
#   h5file, "SFS/fits/1"
# )
# 
