"""
Neural-network version of the ranger random-forest galaxy-property estimator.

Key differences vs. the R/ranger script:
  * Trains on the WHOLE SED (all band magnitudes + redshift) rather than
    engineered nearest-neighbour colors. A NN can learn whatever linear/
    nonlinear combinations of bands it needs internally, so we don't have
    to hand-craft colors the way you do for a tree ensemble. (If you still
    want to give it a head start, set `USE_ALL_PAIRWISE_COLORS = True`
    below to also append every pairwise color, not just adjacent ones.)
  * One multi-output network predicts all 5 targets jointly (shared trunk,
    separate small heads) instead of 5 independent forests. This lets the
    model exploit correlations between StellarMass/SFR10/AGNlum/UV1500/
    UV1500AGN, and is far cheaper to train than 5 separate nets.
  * Standardization of inputs/targets is required for NNs (trees don't
    care about scale, NNs do) -- done with sklearn's StandardScaler,
    fit on the training split only.
  * Early stopping on a held-out validation split (carved out of the
    training set) to avoid overfitting, since NNs don't have the implicit
    regularization of bagging/random feature subsets that RF gets for free.

Requires: torch, numpy, pandas, scikit-learn
"""

import os
import itertools
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
main_dir = os.path.expanduser("~/Documents/Epochs/source_inject/")
fname = "smiles_survey"

depths_path = os.path.expanduser(
    "~/Documents/Epochs/data/depths/smiles_profound_frame_stats.csv"
)
astro_cat_path = os.path.join(main_dir, fname, f"{fname}_mock_training_cat.csv")
phot_cat_path = os.path.join(main_dir, fname, f"{fname}_mock_training_photometry.csv")

BAND_NAMES = [
    "F435W", "F606W", "F775W", "F814W", "F850LP",
    "F090W", "F115W", "F150W", "F200W", "F277W", "F335M", "F356W", "F410M", "F444W",
    "F560W", "F770W", "F1000W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
]
FLUX_COLS = [f"{b}_flux" for b in BAND_NAMES]

TARGET_COLS = ["z", "StellarMass", "SFR10", "AGNlum", "UV1500", "UV1500AGN"]
LOG_TARGET = [False, True, True, True, False, False]

USE_ALL_PAIRWISE_COLORS = False  # set True to also append every color_i_j

SEED = 42
BATCH_SIZE = 512
MAX_EPOCHS = 500
PATIENCE = 20          # early-stopping patience (epochs w/o val improvement)
LR = 1e-3
WEIGHT_DECAY = 1e-5
HIDDEN_DIMS = [256, 128, 64]
DROPOUT = 0.1

DEVICE = torch.device(
    "cuda" if torch.cuda.is_available()
    else "mps" if torch.backends.mps.is_available()
    else "cpu"
)

torch.manual_seed(SEED)
np.random.seed(SEED)

# ----------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------
smiles_depths = pd.read_csv(depths_path)
astro_cat = pd.read_csv(astro_cat_path)
phot_cat = pd.read_csv(phot_cat_path)

train_SED_matrix = phot_cat[FLUX_COLS].to_numpy()  # can have negative fluxes
n_wave = train_SED_matrix.shape[1]

band_noise = smiles_depths[[f"{b}_5sigma_depth_32asec" for b in BAND_NAMES]].to_numpy().ravel()

# asinh ("Luptitude") magnitudes -- same as the R script
mag_matrix = -2.5 * np.log10(np.e) * (
    np.arcsinh(train_SED_matrix / (2 * band_noise)) + np.log(band_noise / 3631.0)
)
mag_cols = [f"band_{i+1}" for i in range(n_wave)]
mag_df = pd.DataFrame(mag_matrix, columns=mag_cols)

# feature_df = pd.concat([astro_cat[["z"]].reset_index(drop=True), mag_df], axis=1)
feature_df = mag_df

if USE_ALL_PAIRWISE_COLORS:
    color_data = {}
    for i, j in itertools.combinations(range(n_wave), 2):
        color_data[f"color_{i+1}_{j+1}"] = mag_matrix[:, i] - mag_matrix[:, j]
    feature_df = pd.concat([feature_df, pd.DataFrame(color_data)], axis=1)

train_parms = astro_cat[TARGET_COLS].copy()
for col, do_log in zip(TARGET_COLS, LOG_TARGET):
    if do_log:
        train_parms[col] = np.log10(train_parms[col])

full_data = pd.concat([feature_df, train_parms], axis=1)

# Drop rows with non-finite values (bad fluxes/targets)
ok = np.isfinite(full_data.to_numpy()).all(axis=1)
n_before = len(full_data)
full_data = full_data.loc[ok].reset_index(drop=True)
print(f"Kept {len(full_data)} / {n_before} rows after removing non-finite cases")

feature_names = list(feature_df.columns)
target_names = TARGET_COLS

X = full_data[feature_names].to_numpy(dtype=np.float32)
Y = full_data[target_names].to_numpy(dtype=np.float32)

# ----------------------------------------------------------------------
# Train / val / test split (80/10/10, same overall train fraction as the R script)
# ----------------------------------------------------------------------
X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=0.2, random_state=SEED
)
X_train, X_val, Y_train, Y_val = train_test_split(
    X_train, Y_train, test_size=0.125, random_state=SEED  # 0.125*0.8 = 0.1 of total
)

x_scaler = StandardScaler().fit(X_train)
y_scaler = StandardScaler().fit(Y_train)

def scale_x(a): return x_scaler.transform(a).astype(np.float32)
def scale_y(a): return y_scaler.transform(a).astype(np.float32)

X_train_s, X_val_s, X_test_s = map(scale_x, (X_train, X_val, X_test))
Y_train_s, Y_val_s = map(scale_y, (Y_train, Y_val))


class SEDDataset(Dataset):
    def __init__(self, X, Y):
        self.X = torch.from_numpy(X)
        self.Y = torch.from_numpy(Y)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]


train_loader = DataLoader(SEDDataset(X_train_s, Y_train_s), batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(SEDDataset(X_val_s, Y_val_s), batch_size=BATCH_SIZE, shuffle=False)

# ----------------------------------------------------------------------
# Model: shared trunk, joint multi-output head
# ----------------------------------------------------------------------
class SEDNet(nn.Module):
    def __init__(self, n_in, n_out, hidden_dims=HIDDEN_DIMS, dropout=DROPOUT):
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


model = SEDNet(n_in=X.shape[1], n_out=len(target_names)).to(DEVICE)
optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.5, patience=8)
loss_fn = nn.MSELoss()

# ----------------------------------------------------------------------
# Train with early stopping
# ----------------------------------------------------------------------
best_val_loss = np.inf
best_state = None
patience_counter = 0

for epoch in range(1, MAX_EPOCHS + 1):
    model.train()
    train_loss = 0.0
    for xb, yb in train_loader:
        xb, yb = xb.to(DEVICE), yb.to(DEVICE)
        optimizer.zero_grad()
        pred = model(xb)
        loss = loss_fn(pred, yb)
        loss.backward()
        optimizer.step()
        train_loss += loss.item() * xb.size(0)
    train_loss /= len(train_loader.dataset)

    model.eval()
    val_loss = 0.0
    with torch.no_grad():
        for xb, yb in val_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            pred = model(xb)
            val_loss += loss_fn(pred, yb).item() * xb.size(0)
    val_loss /= len(val_loader.dataset)
    scheduler.step(val_loss)

    if val_loss < best_val_loss - 1e-5:
        best_val_loss = val_loss
        best_state = {k: v.clone() for k, v in model.state_dict().items()}
        patience_counter = 0
    else:
        patience_counter += 1

    if epoch == 1 or epoch % 10 == 0:
        print(f"epoch {epoch:4d}  train_loss={train_loss:.4f}  val_loss={val_loss:.4f}")

    if patience_counter >= PATIENCE:
        print(f"Early stopping at epoch {epoch} (best val_loss={best_val_loss:.4f})")
        break

model.load_state_dict(best_state)

# ----------------------------------------------------------------------
# Evaluate on held-out test set (in original, unscaled units)
# ----------------------------------------------------------------------
model.eval()
with torch.no_grad():
    pred_test_s = model(torch.from_numpy(X_test_s).to(DEVICE)).cpu().numpy()
pred_test = y_scaler.inverse_transform(pred_test_s)

print("\nTest-set performance:")
for i, targ in enumerate(target_names):
    resid = pred_test[:, i] - Y_test[:, i]
    rmse = np.sqrt(np.mean(resid ** 2))
    ss_res = np.sum(resid ** 2)
    ss_tot = np.sum((Y_test[:, i] - Y_test[:, i].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    print(f"{targ:<12s}  RMSE = {rmse:6.4f}   R^2 = {r2:6.4f}")

# ----------------------------------------------------------------------
# Save model + scalers for later inference
# ----------------------------------------------------------------------
out_dir = os.path.join(main_dir, fname)
os.makedirs(out_dir, exist_ok=True)
torch.save(
    {
        "model_state_dict": model.state_dict(),
        "n_in": X.shape[1],
        "n_out": len(target_names),
        "hidden_dims": HIDDEN_DIMS,
        "dropout": DROPOUT,
        "feature_names": feature_names,
        "target_names": target_names,
        "x_mean": x_scaler.mean_, "x_scale": x_scaler.scale_,
        "y_mean": y_scaler.mean_, "y_scale": y_scaler.scale_,
    },
    os.path.join(out_dir, f"{fname}_sednet.pt"),
)
print(f"\nSaved model + scalers to {os.path.join(out_dir, f'{fname}_sednet.pt')}")