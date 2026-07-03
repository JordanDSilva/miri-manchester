"""
Run inference with a trained SEDNet checkpoint (from train_sed_nn.py) on new
photometry.

Usage (as a script, using the same file layout as training):
    python3 predict_sed_nn.py

Or import the functions and call them on your own photometry DataFrame:
    from predict_sed_nn import load_checkpoint, predict
    ckpt = load_checkpoint("smiles_survey_sednet.pt")
    preds_df = predict(ckpt, phot_df, astro_df["z"].to_numpy(), band_noise)
"""

import os
import itertools
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from astropy.io import fits

# ----------------------------------------------------------------------
# Must match train_sed_nn.py
# ----------------------------------------------------------------------
BAND_NAMES = [
    "F435W", "F606W", "F775W", "F814W", "F850LP",
    "F090W", "F115W", "F150W", "F200W", "F277W", "F335M", "F356W", "F410M", "F444W",
    "F560W", "F770W", "F1000W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
]
FLUX_COLS = [f"{b}_flux" for b in BAND_NAMES]

# which targets were log10'd during training -- keyed by name so order-safe
LOG_TARGET_MAP = {
    "z" : False,
    "StellarMass": True,
    "SFR10": True,
    "AGNlum": True,
    "UV1500": False,
    "UV1500AGN": False,
}

DEVICE = torch.device(
    "cuda" if torch.cuda.is_available()
    else "mps" if torch.backends.mps.is_available()
    else "cpu"
)


class SEDNet(nn.Module):
    def __init__(self, n_in, n_out, hidden_dims, dropout):
        super().__init__()
        layers = []
        prev = n_in
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        self.trunk = nn.Sequential(*layers)
        self.head = nn.Linear(prev, n_out)

    def forward(self, x):
        return self.head(self.trunk(x))


def load_checkpoint(pt_path):
    ckpt = torch.load(pt_path, map_location=DEVICE, weights_only=False)
    model = SEDNet(
        n_in=ckpt["n_in"],
        n_out=ckpt["n_out"],
        hidden_dims=ckpt["hidden_dims"],
        dropout=ckpt["dropout"],
    ).to(DEVICE)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    ckpt["model"] = model
    return ckpt


def _build_feature_df(flux_matrix, redshift, band_noise, feature_names):
    """Recreate the same features used at training time (mags + optionally
    every pairwise color), then select/order columns to match feature_names
    exactly -- so this works regardless of whether USE_ALL_PAIRWISE_COLORS
    was on or off during training."""
    n_wave = flux_matrix.shape[1]

    mag_matrix = -2.5 * np.log10(np.e) * (
        np.arcsinh(flux_matrix / (2 * band_noise)) + np.log(band_noise / 3631.0)
    )
    mag_cols = [f"band_{i+1}" for i in range(n_wave)]
    df = pd.DataFrame(mag_matrix, columns=mag_cols)
    df.insert(0, "z", redshift)

    needs_colors = any(c.startswith("color_") for c in feature_names)
    if needs_colors:
        color_data = {}
        for i, j in itertools.combinations(range(n_wave), 2):
            name = f"color_{i+1}_{j+1}"
            if name in feature_names:
                color_data[name] = mag_matrix[:, i] - mag_matrix[:, j]
        df = pd.concat([df, pd.DataFrame(color_data)], axis=1)

    missing = [c for c in feature_names if c not in df.columns]
    if missing:
        raise ValueError(f"Could not reconstruct required features: {missing}")

    return df[feature_names].to_numpy(dtype=np.float32)

def _build_feature_dfV2(flux_matrix, band_noise, feature_names):
    """Recreate the same features used at training time (mags + optionally
    every pairwise color), then select/order columns to match feature_names
    exactly -- so this works regardless of whether USE_ALL_PAIRWISE_COLORS
    was on or off during training."""
    n_wave = flux_matrix.shape[1]

    mag_matrix = -2.5 * np.log10(np.e) * (
        np.arcsinh(flux_matrix / (2 * band_noise)) + np.log(band_noise / 3631.0)
    )
    mag_cols = [f"band_{i+1}" for i in range(n_wave)]
    df = pd.DataFrame(mag_matrix, columns=mag_cols)

    needs_colors = any(c.startswith("color_") for c in feature_names)
    if needs_colors:
        color_data = {}
        for i, j in itertools.combinations(range(n_wave), 2):
            name = f"color_{i+1}_{j+1}"
            if name in feature_names:
                color_data[name] = mag_matrix[:, i] - mag_matrix[:, j]
        df = pd.concat([df, pd.DataFrame(color_data)], axis=1)

    missing = [c for c in feature_names if c not in df.columns]
    if missing:
        raise ValueError(f"Could not reconstruct required features: {missing}")

    return df[feature_names].to_numpy(dtype=np.float32)


def predict(ckpt, phot_df, redshift, band_noise, flux_cols=FLUX_COLS):
    """
    ckpt        : output of load_checkpoint()
    phot_df     : DataFrame containing flux_cols (e.g. 'F435W_flux', ...)
    redshift    : array-like of length nrow(phot_df) -- known/estimated z
    band_noise  : array-like, 5-sigma depths per band, same order as BAND_NAMES

    Returns a DataFrame of predicted target values in original physical
    units (i.e. StellarMass/SFR10/AGNlum already un-logged back from log10).
    """
    flux_matrix = phot_df[flux_cols].to_numpy(dtype=np.float64)
    band_noise = np.asarray(band_noise, dtype=np.float64)
    redshift = np.asarray(redshift, dtype=np.float64)

    # X = _build_feature_df(flux_matrix, redshift, band_noise, ckpt["feature_names"])
    X = _build_feature_dfV2(flux_matrix, band_noise, ckpt["feature_names"])


    finite_mask = np.isfinite(X).all(axis=1)
    if not finite_mask.all():
        print(f"Warning: {(~finite_mask).sum()} rows have non-finite features and will be NaN in output")

    X_s = (X - ckpt["x_mean"]) / ckpt["x_scale"]
    X_s = np.nan_to_num(X_s, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)

    with torch.no_grad():
        pred_s = ckpt["model"](torch.from_numpy(X_s).to(DEVICE)).cpu().numpy()

    pred = pred_s * ckpt["y_scale"] + ckpt["y_mean"]

    out = pd.DataFrame(pred, columns=ckpt["target_names"])
    for col in out.columns:
        if LOG_TARGET_MAP.get(col, False):
            out[col] = 10 ** out[col]

    out.loc[~finite_mask, :] = np.nan
    return out


# ----------------------------------------------------------------------
# Example usage matching the same file layout as train_sed_nn.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main_dir = os.path.expanduser("~/Documents/Epochs/source_inject/")
    fname = "smiles_survey"
    sednet_fname = "smiles_survey"

    depths_path = os.path.expanduser(
        "~/Documents/Epochs/data/depths/smiles_profound_frame_stats.csv"
    )
    astro_cat_path = os.path.join(main_dir, fname, f"{fname}_mock_cat.csv")
    phot_cat_path = os.path.join(main_dir, fname, f"{fname}_mock_eazy_photometry.csv")
    pt_path = os.path.join(main_dir, sednet_fname, f"{sednet_fname}_sednet.pt")

    with fits.open("/Users/jordan/Documents/Epochs/source_inject/smiles_survey/eazy/smiles_survey_zphot_agn_blue_sfhz.out.zout.fits") as hdul:
        data = hdul[1].data
        redshifts = np.array(data["z_ml"])

    smiles_depths = pd.read_csv(depths_path)
    astro_cat = pd.read_csv(astro_cat_path)
    phot_cat = pd.read_csv(phot_cat_path)
    band_noise = smiles_depths[[f"{b}_5sigma_depth_32asec" for b in BAND_NAMES]].to_numpy().ravel()

    ckpt = load_checkpoint(pt_path)
    # preds = predict(ckpt, phot_cat, astro_cat["z"].to_numpy(), band_noise)
    preds = predict(ckpt, phot_cat, redshifts, band_noise)

    print(preds.head())
    preds.to_csv(os.path.join(main_dir, fname, f"{fname}_nn_predictions.csv"), index=False)
    print(f"Saved predictions to {os.path.join(main_dir, fname, f'{fname}_nn_predictions.csv')}")