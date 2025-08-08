#!/usr/bin/env python3
"""
reconstruction.py

Utilities for simulation, calibration (LUT building), reconstruction and visualization
for a 2D intensity-ratio detector (5x5 CsI tiles read by 4 SiPM boards).

Usage examples (from repo root):

# 1) Simulate a beam scan and build calibration LUT
python software/analysis/reconstruction.py --simulate_scan --scan_positions 6 --out lut.npy

# 2) Use LUT to reconstruct measurements saved in CSV
python software/analysis/reconstruction.py --reconstruct --lut lut.npy --csv data/daq_sample.csv --out reconstructed.npy

# 3) Plot reconstruction heatmap from reconstructed positions
python software/analysis/reconstruction.py --plot_recon --recon reconstructed.npy

Author: Nitesh Kumar & Abhishek Keshri (NIT Patna)
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.stats import binned_statistic_2d
import pandas as pd
import joblib
import json

# -------------------------
# Configuration / Defaults
# -------------------------
DETECTOR_SIZE_CM = 5.0  # physical size of 5x5 grid (example): 5 cm -> tiles are 1 cm
GRID_N = 5              # 5x5 tiles
NOISE_STD = 0.02        # relative noise
ATTENUATION_FACTOR = 0.8  # how much signals attenuate across tiles
RANDOM_SEED = 42


def safe_mkdir(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def ratios_from_signals(E1, E2, E3, E4):
    """
    Compute the two ratio features used in the proposal:
      X_ratio = E2 / (E2 + E3)   (vertical pair)
      Y_ratio = E1 / (E1 + E4)   (horizontal pair)
    Returns (x_ratio, y_ratio)
    Handles division by zero by returning NaN.
    """
    denom_v = (E2 + E3)
    denom_h = (E1 + E4)
    xr = np.where(denom_v != 0, E2 / denom_v, np.nan)
    yr = np.where(denom_h != 0, E1 / denom_h, np.nan)
    return xr, yr



# Simulator

class DetectorSimulator:
    """
    Simulate SiPM signals (E1..E4) for hits at given physical (x,y) positions (in cm).
    The position coordinate system: origin at detector center, x to right, y up.
    """

    def __init__(self, size_cm=DETECTOR_SIZE_CM, grid_n=GRID_N,
                 base_signal=1.0, noise_std=NOISE_STD, attenuation=ATTENUATION_FACTOR, seed=RANDOM_SEED):
        self.size_cm = size_cm
        self.grid_n = grid_n
        self.base_signal = base_signal
        self.noise_std = noise_std
        self.attenuation = attenuation
        self.rng = np.random.default_rng(seed)

    def tile_centers(self):
        """Return centers of tiles in physical coordinates (x,y) arrays, origin at center."""
        tile_size = self.size_cm / self.grid_n
        coords_1d = (np.arange(self.grid_n) - (self.grid_n - 1) / 2.0) * tile_size
        xs, ys = np.meshgrid(coords_1d, coords_1d)
        return xs.flatten(), ys.flatten()

    def simulate_hit(self, x_cm, y_cm, spread=0.4):
        """
        Given a hit at (x_cm, y_cm), simulate four SiPM integrated signals E1..E4.
        Model: each SiPM board receives light that decays with distance from hit
        and also has cross-coupling from neighbor tiles. This is a simplified model.
        Returns E1, E2, E3, E4 as floats.
        Board orientations:
          - E1: top (north)
          - E2: right (east)
          - E3: left (west)
          - E4: bottom (south)
        """
        # compute normalized distance of hit to each board center (approx)
        # For simplicity model: board signal ~ base * exp(-d^2/(2*sigma^2)), but sigma tuned per board
        # Place board centers outside the tile area at midpoints of sides
        half = self.size_cm / 2.0
        # Board centers
        cE1 = np.array([0.0, half + 0.5])   # top
        cE2 = np.array([half + 0.5, 0.0])   # right
        cE3 = np.array([-half - 0.5, 0.0])  # left
        cE4 = np.array([0.0, -half - 0.5])  # bottom

        p = np.array([x_cm, y_cm])
        def board_signal(center):
            d2 = np.sum((p - center) ** 2)
            # baseline Gaussian coupling + small attenuation factor
            s = self.base_signal * np.exp(-d2 / (2 * spread * spread))
            return s

        s1 = board_signal(cE1)
        s2 = board_signal(cE2)
        s3 = board_signal(cE3)
        s4 = board_signal(cE4)

        # apply inter-board attenuation / bias to mimic geometry
        s1 *= 1.0
        s2 *= 1.0
        s3 *= 1.0
        s4 *= 1.0

        # Add multiplicative and additive noise
        noise = self.rng.normal(1.0, self.noise_std, size=4)
        signals = np.array([s1, s2, s3, s4]) * noise

        # Add small baseline / dark counts
        baseline = self.rng.normal(0.01, 0.005, size=4)
        signals = signals + baseline
        return signals.tolist()

    def simulate_beam_scan(self, n_positions=6, samples_per_pos=500, out_csv=None):
        """
        Simulate a beam scan across x positions (center row) and produce CSV-like pandas DataFrame:
        columns: ['E1','E2','E3','E4','x_true','y_true','timestamp']
        """
        xs = np.linspace(-self.size_cm/2 + 0.5, self.size_cm/2 - 0.5, n_positions)
        y_true = 0.0  # scanning along x only for beam test
        rows = []
        t = 0.0
        for x in xs:
            for i in range(samples_per_pos):
                E1, E2, E3, E4 = self.simulate_hit(x, y_true)
                rows.append([E1, E2, E3, E4, x, y_true, t])
                t += 0.001
        df = pd.DataFrame(rows, columns=['E1','E2','E3','E4','x_true','y_true','timestamp'])
        if out_csv:
            safe_mkdir(out_csv)
            df.to_csv(out_csv, index=False)
            print(f"Saved simulated beam scan to {out_csv}")
        return df



class LUTCalibrator:
    """
    Build a 2D interpolant (LUT) mapping (x_ratio, y_ratio) -> (x_cm, y_cm).
    Uses sampled calibration points (from beam or controlled illumination).
    """

    def __init__(self, method='linear'):
        self.method = method
        self.interpolant_x = None
        self.interpolant_y = None

    def build_from_dataframe(self, df):
        """
        df must contain columns E1,E2,E3,E4,x_true,y_true
        """
        E1 = df['E1'].values
        E2 = df['E2'].values
        E3 = df['E3'].values
        E4 = df['E4'].values
        xr, yr = ratios_from_signals(E1, E2, E3, E4)
        valid = ~(np.isnan(xr) | np.isnan(yr))
        pts = np.vstack([xr[valid], yr[valid]]).T
        x_targets = df['x_true'].values[valid]
        y_targets = df['y_true'].values[valid]

        # Build interpolants
        if self.method == 'linear':
            self.interpolant_x = LinearNDInterpolator(pts, x_targets, fill_value=np.nan)
            self.interpolant_y = LinearNDInterpolator(pts, y_targets, fill_value=np.nan)
        else:
            # fallback to nearest if linear fails for some points
            self.interpolant_x = NearestNDInterpolator(pts, x_targets)
            self.interpolant_y = NearestNDInterpolator(pts, y_targets)

        print("LUT built with method:", self.method)
        return self

    def save(self, path):
        safe_mkdir(path)
        joblib.dump({'method': self.method,
                     'interpolant_x': self.interpolant_x,
                     'interpolant_y': self.interpolant_y}, path)
        print("Saved LUT to", path)

    @staticmethod
    def load(path):
        data = joblib.load(path)
        obj = LUTCalibrator(method=data.get('method', 'linear'))
        obj.interpolant_x = data['interpolant_x']
        obj.interpolant_y = data['interpolant_y']
        print("Loaded LUT from", path)
        return obj

    def lookup(self, E1, E2, E3, E4):
        xr, yr = ratios_from_signals(np.array([E1]), np.array([E2]), np.array([E3]), np.array([E4]))
        if np.isnan(xr[0]) or np.isnan(yr[0]):
            return (np.nan, np.nan)
        x_pred = float(self.interpolant_x(xr[0], yr[0]))
        y_pred = float(self.interpolant_y(xr[0], yr[0]))
        return (x_pred, y_pred)



def reconstruct_dataframe_with_lut(df, lut: LUTCalibrator):
    E1 = df['E1'].values
    E2 = df['E2'].values
    E3 = df['E3'].values
    E4 = df['E4'].values
    xr, yr = ratios_from_signals(E1, E2, E3, E4)
    preds_x = lut.interpolant_x(xr, yr)
    preds_y = lut.interpolant_y(xr, yr)
    out = df.copy()
    out['xr'] = xr
    out['yr'] = yr
    out['x_rec'] = preds_x
    out['y_rec'] = preds_y
    return out



def parse_daq_csv(path, require_columns=None):
    """
    Parse a CSV produced by DAQ or serial logger.
    Expect columns: E1,E2,E3,E4,... optionally x_true,y_true
    """
    df = pd.read_csv(path)
    if require_columns:
        missing = [c for c in require_columns if c not in df.columns]
        if missing:
            raise ValueError(f"Missing columns in DAQ CSV: {missing}")
    return df



# Visualization
def plot_2d_histogram(df, value_col='x_rec', bins=100, title='Reconstructed positions'):
    """
    Simple 2D histogram plot of reconstructed (x_rec,y_rec)
    """
    x = df['x_rec'].values
    y = df['y_rec'].values
    mask = ~(np.isnan(x) | np.isnan(y))
    if np.sum(mask) == 0:
        print("No valid reconstructed points to plot.")
        return
    hb = plt.hist2d(x[mask], y[mask], bins=bins)
    plt.colorbar()
    plt.xlabel('x (cm)')
    plt.ylabel('y (cm)')
    plt.title(title)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def plot_ratio_histograms(df):
    xr = df['xr'].values
    yr = df['yr'].values
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.hist(xr[~np.isnan(xr)], bins=50)
    plt.title('X ratio distribution')
    plt.subplot(1,2,2)
    plt.hist(yr[~np.isnan(yr)], bins=50)
    plt.title('Y ratio distribution')
    plt.show()



def main_cli():
    parser = argparse.ArgumentParser(description="Calibration & Reconstruction utilities for 2D SiPM detector")
    parser.add_argument("--simulate_scan", action="store_true", help="Simulate beam scan dataset")
    parser.add_argument("--scan_positions", type=int, default=6, help="Number of beam positions for simulated scan")
    parser.add_argument("--samples_per_pos", type=int, default=500, help="Samples per simulated position")
    parser.add_argument("--out", type=str, default="data/simulated_scan.csv", help="Output CSV or LUT file")
    parser.add_argument("--build_lut", action="store_true", help="Build LUT from CSV file (requires --out to be lut path)")
    parser.add_argument("--csv", type=str, default=None, help="CSV file with E1..E4 columns for reconstruction")
    parser.add_argument("--lut", type=str, default="data/lut.pkl", help="LUT file path")
    parser.add_argument("--reconstruct", action="store_true", help="Reconstruct CSV using LUT")
    parser.add_argument("--plot_recon", action="store_true", help="Plot reconstructed results from NPY/CSV")
    args = parser.parse_args()

    if args.simulate_scan:
        sim = DetectorSimulator()
        df = sim.simulate_beam_scan(n_positions=args.scan_positions, samples_per_pos=args.samples_per_pos, out_csv=args.out)
        print(df.head())
        return

    if args.build_lut:
        if not args.csv:
            raise SystemExit("Please provide --csv with calibration CSV (E1..E4,x_true,y_true)")
        df_cal = parse_daq_csv(args.csv, require_columns=['E1','E2','E3','E4','x_true','y_true'])
        lut = LUTCalibrator(method='linear').build_from_dataframe(df_cal)
        lut.save(args.out)
        return

    if args.reconstruct:
        if not args.csv:
            raise SystemExit("Please provide --csv with data to reconstruct (E1..E4)")
        lut = LUTCalibrator.load(args.lut)
        df = parse_daq_csv(args.csv, require_columns=['E1','E2','E3','E4'])
        out = reconstruct_dataframe_with_lut(df, lut)
        # Save results
        safe_mkdir(args.out)
        out.to_csv(args.out, index=False)
        print(f"Saved reconstructed CSV to {args.out}")
        return

    if args.plot_recon:
        if not args.csv:
            raise SystemExit("Please provide --csv with reconstructed CSV (must have x_rec,y_rec)")
        df = parse_daq_csv(args.csv)
        plot_2d_histogram(df)
        return

    parser.print_help()


if __name__ == "__main__":
    main_cli()
