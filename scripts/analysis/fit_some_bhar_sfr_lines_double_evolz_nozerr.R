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
library(dplyr)

set.seed(666)

prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

BHARSFR_df = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/smbh/prospect_robust_BHAR_SFR_ratio.csv"))

LAGN2BHAR = function(epsilon){( 3.154e7 * 1e-7 * (epsilon * 3e8^2)**-1 / 2e30)}
epsilon_samples = runif(n = 100000, min = 0.06, max = 0.4)

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

# --- Data ---
z_obs = prospect_sample$zQ50
x_obs = log10(prospect_sample$StellarMassQ50)
y_obs = log10(BHARSFR_df_sample$BHARSFR10Q50)

zerr = (prospect_sample$zQ84 - prospect_sample$zQ16) * 0.5
xerr = (log10(prospect_sample$StellarMassQ84) - log10(prospect_sample$StellarMassQ16)) * 0.5
yerr = (log10(BHARSFR_df_sample$BHARSFR10Q84) - log10(BHARSFR_df_sample$BHARSFR10Q16)) * 0.5


cutoff_mstar = 0.0
fit_idx = x_obs >= cutoff_mstar & prospect_sample$FitType == "StellarAGN"

N = length(which(fit_idx))

covarray = array(0, dim = c(2, 2, N))
for (j in 1:N){
  read_stub = paste0(
    "~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/",
    "/galID_", prospect_sample$CATID[which(fit_idx)[j]], ".rds"
  )
  read_data = readRDS(read_stub)
  post = read_data$post[, c("StellarMass", "SFR10", "AGNlum", "z")]
  
  Npost = nrow(post)
  
  x = as.numeric(log10(post[, "StellarMass"]))
  y = as.numeric(log10(post[, "SFR10"] / ( LAGN2BHAR(epsilon_samples[sample(1:Npost)]) * post[, "AGNlum"])))
  
  ok = is.finite(x) & is.finite(y)
  
  if (sum(ok) < 2){
    # not enough valid samples to compute covariance — fall back to something safe
    warning(sprintf("j = %d: insufficient finite samples (%d), using diagonal fallback", j, sum(ok)))
    covarray[, , j] = diag(c(
      var(x[is.finite(x)], na.rm = TRUE),
      var(y[is.finite(y)], na.rm = TRUE)
    ))
    next
  }
  
  covarray[, , j] = cov(cbind(x[ok], y[ok]))
}

Data = list()
Data$z_obs = z_obs[fit_idx]; Data$x_obs = x_obs[fit_idx]; Data$y_obs = y_obs[fit_idx]
Data$N = N; Data$covarray = covarray

# meat 
z_piv = 3

xmstar = seq(-15, 15, 0.01)
dbl_pl_vec = function(x, M0, alpha, beta, S0){
  xx = 10^(x - M0)
  S0 - log10(xx^(-alpha) + xx^(-beta))
}

# vectorised in z; returns an N x 5 matrix
redshift_evol_parms = function(z, parm){
  zc = z - z_piv
  quad = function(i) parm[i] + parm[i+1]*zc + parm[i+2]*zc^2
  cbind(M0 = quad(1), alpha = quad(4), logdelta = quad(7),
        S0 = quad(10), logsint = quad(13))
}
redshift_evol_parms_for_dbl = function(z, parm){
  zc = z - z_piv
  quad = function(i) parm[i] + parm[i+1]*zc + parm[i+2]*zc^2
  cbind(M0 = quad(1), alpha = quad(4), beta = quad(4) - exp(quad(7)),
        S0 = quad(10), logsint = quad(13))
}

Model = function(parm, Data){
  N     = Data$N
  coef  = parm[1:15]
  xtrue = parm[16:(15 + N)]
  
  P = redshift_evol_parms(Data$z_obs, coef)
  M0 = P[,1]; alpha = P[,2]; delta = exp(P[,3]); S0 = P[,4]
  sint2 = exp(2 * P[,5])
  beta  = alpha - delta
  
  pm = c(7,0,0,   1,0,0,   0,0,0,   2,0,0,   log(0.3),0,0)
  ps = c(3,1,0.5, 2,1,0.5, 1,0.5,0.25, 3,1,0.5, 1,0.5,0.25)
  # LPrior = sum(dnorm(coef, pm, ps, log = TRUE))
  LPrior = 0
  
  y_model = dbl_pl_vec(xtrue, M0, alpha, beta, S0)
  
  eps = 1e-3
  S   = Data$covarray                     # 3 x 3 x N, order (z, x, y)
  sx2 = S[1,1,] + eps
  sy2 = S[2,2,] + eps + sint2             # intrinsic scatter in y only
  sxy = S[1,2,]
  
  dx = Data$x_obs - xtrue
  dy = Data$y_obs - y_model
  
  detS = sx2*sy2 - sxy^2
  if (any(!is.finite(detS)) || any(detS <= 0)) {
    return(list(LP = -1e100, Dev = 1e100, Monitor = c(-1e100, coef),
                yhat = rep(NA_real_, N), parm = parm))
  }
  
  quad = (dx^2*sy2 - 2*dx*dy*sxy + dy^2*sx2) / detS
  LL   = sum(-log(2*pi) - 0.5*log(detS) - 0.5*quad)
  
  LP = LL + LPrior
  
  list(LP = LP, Dev = -2*LL, Monitor = c(LP, coef), yhat = y_model, parm = parm)
}

nm = c("M0","alpha","logdelta","S0","logsint")
Data$parm.names = c(paste0(rep(nm, each = 3), "_c", 1:3),
                    paste0("xtrue", 1:N))
Data$mon.names  = c("LP", paste0(rep(nm, each = 3), "_c", 1:3))

PGF = function(Data){
  c(7,0,0,  1,0,0,  0,0,0,  2,0,0,  log(0.3),0,0,
    # Data$z_obs + rnorm(Data$N, 0, 0.02),
    Data$x_obs + rnorm(Data$N, 0, 0.05))
}

Initial.Values = PGF(Data)

Fit = LaplacesDemon(
  Model = Model,
  Data  = Data,
  Initial.Values = Initial.Values,
  Covar = NULL,
  Iterations = 100000,
  Status = 1000,
  Thinning = 10,
  Algorithm = "CHARM",       # adaptive algorithm, good default for this kind of problem
  Specs = NULL
)

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
zvec = seq(1, 10, 0.1)
parm_z_func = function(z, parm){
  zc = z - z_piv
  parm[1] + parm[2]*zc + parm[3]*zc^2
} 

post <- Fit$Monitor[, paste0(rep(nm, each = 3), "_c", 1:3)]
df_quantiles = colQuantiles(
  x = post,
  probs = c(0.5, 0.16, 0.84),
  na.rm = TRUE
)
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

zgrid = seq(1, 10, 0.1)
M0z_samples = t(apply(as.matrix(post[, c("M0_c1", "M0_c2", "M0_c3")]), 1, function(p) parm_z_func(zgrid, p)))
M0z_q = t(apply(M0z_samples, 2, quantile, probs = c(0.50, 0.16, 0.84)))

alphaz_samples = t(apply(as.matrix(post[, c("alpha_c1", "alpha_c2", "alpha_c3")]), 1, function(p) parm_z_func(zgrid, p)))
alphaz_q = t(apply(alphaz_samples, 2, quantile, probs = c(0.50, 0.16, 0.84)))

logdeltaz_samples = t(apply(as.matrix(post[, c("logdelta_c1", "logdelta_c2", "logdelta_c3")]), 1, function(p) parm_z_func(zgrid, p)))
logdeltaz_q = t(apply(logdeltaz_samples, 2, quantile, probs = c(0.50, 0.16, 0.84)))

betaz_samples = alphaz_samples - exp(logdeltaz_samples)
betaz_q = t(apply(betaz_samples, 2, quantile, probs = c(0.50, 0.16, 0.84)))

S0z_samples = t(apply(as.matrix(post[, c("S0_c1", "S0_c2", "S0_c3")]), 1, function(p) parm_z_func(zgrid, p)))
S0z_q = t(apply(S0z_samples, 2, quantile, probs = c(0.50, 0.16, 0.84)))

sint_samples = exp(t(apply(as.matrix(post[, c("logsint_c1", "logsint_c2", "logsint_c3")]), 1, function(p) parm_z_func(zgrid, p))))
sint_q = t(apply(sint_samples, 2, quantile, probs = c(0.50, 0.16, 0.84)))

magplot(
  zgrid,
  M0z_q[,1],
  ylim = c(4,12)
)
lines(
  zgrid,
  M0z_q[,2]
)
lines(
  zgrid,
  M0z_q[,3]
)

magplot(
  zgrid,
  alphaz_q[,1],
  # ylim = c(4,12)
)
lines(
  zgrid,
  alphaz_q[,2]
)
lines(
  zgrid,
  alphaz_q[,3]
)

magplot(
  zgrid,
  betaz_q[,1],
  ylim = c(-0.5,1),
  xlim = c(1,10)
)
lines(
  zgrid,
  betaz_q[,2]
)
lines(
  zgrid,
  betaz_q[,3]
)

magplot(
  zgrid,
  S0z_q[,1],
)
lines(
  zgrid,
  S0z_q[,2]
)
lines(
  zgrid,
  S0z_q[,3]
)

df_fit = list()
for(i in 1:length(zmids)){
  
  zidx = z_obs >= zbins[i] & z_obs < zbins[i + 1] & fit_idx
  # --- Extract posterior on model parameters only (drop latent x_true columns) ---
  
  # Posterior samples of the model parameters
  
  sint_z = apply(post, 1, function(post_row) {
    dbl_test_parm <- redshift_evol_parms_for_dbl(zmids[i], post_row)
    exp(dbl_test_parm[5])
  })
  
  y_samples <- apply(post, 1, function(post_row) {
    dbl_test_parm <- redshift_evol_parms_for_dbl(zmids[i], post_row)
    dbl_pl_vec(
      x_grid,
      dbl_test_parm[1],
      dbl_test_parm[2],
      dbl_test_parm[3],
      dbl_test_parm[4]
    )
  })
  
  # Add intrinsic scatter to each posterior realisation
  y_pred_samples <- y_samples
  
  for (j in seq_len(nrow(post))) {
    y_pred_samples[, j] <- y_samples[, j] +
      rnorm(length(x_grid), mean = 0, sd = (sint_z[j]))
  }
  
  # Quantiles of the fitted mean
  y_quantiles <- t(apply(
    y_samples, 1, quantile,
    probs = c(0.50, 0.16, 0.84),
    na.rm = TRUE
  ))
  
  # Quantiles of the predictive distribution
  y_pred_quantiles <- t(apply(
    y_pred_samples, 1, quantile,
    probs = c(0.50, 0.16, 0.84),
    na.rm = TRUE
  ))
  
  colnames(y_quantiles) <- c("Q50", "Q16", "Q84")
  colnames(y_pred_quantiles) <- c("Q50", "Q16", "Q84")
  
  # # --- Apply across all posterior samples ---
  Mbreak_quantiles = apply(y_pred_quantiles, 2, function(x){find_zero_from_grid(x_grid, y_vals = x)})
  Mbreak_quantiles[c("Q16", "Q84")] = Mbreak_quantiles[c("Q84", "Q16")]
  Mbreak_quantiles[is.na(Mbreak_quantiles)] = -999
  Mbreak_quantiles["lo"] = Mbreak_quantiles["Q50"] - Mbreak_quantiles["Q16"]
  Mbreak_quantiles["hi"] = Mbreak_quantiles["Q84"] - Mbreak_quantiles["Q50"]
  
  y_quantiles = cbind(y_quantiles, "x" = x_grid)
  y_pred_quantiles = cbind(y_pred_quantiles, "x" = x_grid)
  
  # --- Base plot ---
  suppressWarnings(dir.create(
    paste0(
      "~/Documents/JWST-MIRI-AGN/plots/rough-plots/fit_bhar_sfr_fitz", "/", cutoff_mstar, "/"
    ), recursive = TRUE
  ))
  pdf(
    paste0(
      "~/Documents/JWST-MIRI-AGN/plots/rough-plots/fit_bhar_sfr_fitz", "/", cutoff_mstar, "/z", i, ".pdf"
    ),
    width = 8,
    height = 6
  )
  magplot(
    (x_obs[zidx]),
    (y_obs[zidx]),
    xlim = c(1, 12),
    ylim = c(-5, 5),
    xlab = "log10(Stellar Mass)",
    ylab = "log10(BHAR/SFR)"
  )
  magerr(
    (x_obs[zidx]),
    (y_obs[zidx]),
    xlo = xerr[zidx],
    ylo = yerr[zidx]
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
  polygon(
    x = c(x_grid, rev(x_grid)),
    y = c(y_pred_quantiles[, "Q16"], rev(y_pred_quantiles[, "Q84"])),
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
  
  # # --- Mbreak ---
  points(
    Mbreak_quantiles["Q50"],
    0,
    pch = 17,
    col = "red",
    cex = 1.5
  )
  dev.off()
  
  df_fit[[i]] = list(
    "fit" = y_quantiles,
    "pred_fit" = y_pred_quantiles,
    "Mbreak" = Mbreak_quantiles
  )
}

df_fit = df_fit[!sapply(df_fit, is.null)]
# saveRDS(
#   df_fit,
#   "~/Documents/JWST-MIRI-AGN/data/prospect/prospect_robust_SMBH_evolz_results.rds"
# )
# df_fit = readRDS("~/Documents/JWST-MIRI-AGN/data/prospect/prospect_robust_SMBH_evolz_results.rds")


parm_zevol_df = cbind(
  "z" = zvec,
  M0z_quantiles,
  alphaz_quantiles,
  betaz_quantiles,
  S0z_quantiles,
  sint2z_quantiles
)


h5file <- path.expand("~/Documents/JWST-MIRI-AGN/data/analysis_data.h5")
# Create the file if it doesn't exist yet; otherwise we append into it
if (!file.exists(h5file)) {
  h5createFile(h5file)
}

h5ls(h5file)
if (any(h5ls(h5file)$name == "SFRBHAR")) {
  h5delete(h5file, "SFRBHAR")
}
# h5delete(h5file, "SFRBHAR")
h5createGroup(h5file, "SFRBHAR")
h5write(df_cat, file = h5file, name = paste0("SFRBHAR/Params"))
h5write(df_cat, file = h5file, name = paste0("SFRBHAR/Mbreak"))

h5createGroup(h5file, "SFRBHAR/fits")
for(i in 1:length(df_fit)){
  h5write(
    data.frame(df_fit[[i]]$fit),
    file = h5file,
    name = paste0("SFRBHAR/fits/", i)
  )
}

h5createGroup(h5file, "SFRBHAR/pred_fits")
for(i in 1:length(df_fit)){
  h5write(
    data.frame(df_fit[[i]]$pred_fit),
    file = h5file,
    name = paste0("SFRBHAR/pred_fits/", i)
  )
}


parm_z_func = function(z, parm){
  zc = z - z_piv
  parm[1] + parm[2]*zc + parm[3]*zc^2
} 
magplot(
  zvec,
  parm_z_func(zvec, as.numeric(df_cat[c("M0_c1Q50", "M0_c2Q50", "M0_c3Q50")])),
  type = "l"
)
magplot(
  zvec,
  parm_z_func(zvec, as.numeric(df_cat[c("alpha_c1Q50", "alpha_c2Q50", "alpha_c3Q50")])),
  type = "l"
)
magplot(
  zvec,
  parm_z_func(zvec, as.numeric(df_cat[c("beta_c1Q50", "beta_c2Q50", "beta_c3Q50")])),
  type = "l"
)
magplot(
  zvec,
  parm_z_func(zvec, as.numeric(df_cat[c("S0_c1Q50", "S0_c2Q50", "S0_c3Q50")])),
  type = "l"
)
# lines(
#   zvec,
#   parm_z_func(zvec, as.numeric(df_cat[c("M0_c1Q16", "M0_c2Q16", "M0_c3Q16")])),
# )
