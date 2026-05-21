"""Reproduce the R-alpha sensitivity analysis figure."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.dual_constraint_idw import (  # noqa: E402
    IDWConfig,
    coordinates_from_df,
    identify_anomalies,
    improved_idw_predict,
    load_station_data,
    pairwise_distances,
    values_from_df,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default=ROOT / "data" / "example_rainfall_stations.csv")
    parser.add_argument("--out-dir", default=ROOT / "outputs")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_station_data(args.data)
    coords = coordinates_from_df(df)
    values = values_from_df(df)
    dist_matrix = pairwise_distances(coords)
    base_config = IDWConfig()
    anomaly_flags, _ = identify_anomalies(
        values,
        dist_matrix,
        n_neighbors=base_config.anomaly_neighbors,
        threshold=base_config.anomaly_threshold,
    )

    r_values = [20, 30, 40, 50]
    alpha_values = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
    rows = []

    for radius in r_values:
        for alpha in alpha_values:
            config = IDWConfig(radius_km=float(radius), alpha=float(alpha))
            errors = []
            for i in range(len(values)):
                train_idx = np.delete(np.arange(len(values)), i)
                pred = improved_idw_predict(
                    coords[train_idx],
                    values[train_idx],
                    anomaly_flags[train_idx],
                    coords[i],
                    config,
                )
                errors.append(pred - values[i])
            errors = np.asarray(errors)
            rows.append({
                "R_km": radius,
                "alpha": alpha,
                "RMSE_mm_h": float(np.sqrt(np.mean(errors**2))),
                "MAE_mm_h": float(np.mean(np.abs(errors))),
            })

    results = pd.DataFrame(rows)
    results.to_csv(out_dir / "fig5_sensitivity_results.csv", index=False)

    width, height = 1800, 1350
    margin_l, margin_r, margin_t, margin_b = 220, 120, 90, 190
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    colors = {20: (205, 102, 29), 30: (34, 139, 34), 40: (178, 34, 34), 50: (104, 34, 139)}

    x_min, x_max = 0.0, 1.0
    y_min = float(results["RMSE_mm_h"].min()) - 0.08
    y_max = float(results["RMSE_mm_h"].max()) + 0.08

    def xp(x):
        return int(margin_l + (x - x_min) / (x_max - x_min) * plot_w)

    def yp(y):
        return int(margin_t + (y_max - y) / (y_max - y_min) * plot_h)

    # Grid and axes.
    for x_tick in np.linspace(0, 1, 6):
        xpix = xp(float(x_tick))
        draw.line([(xpix, margin_t), (xpix, margin_t + plot_h)], fill=(220, 220, 220), width=1)
        draw.text((xpix - 14, margin_t + plot_h + 18), f"{x_tick:.1f}", fill="black", font=font)
    for y_tick in np.linspace(round(y_min, 1), round(y_max, 1), 6):
        ypix = yp(float(y_tick))
        draw.line([(margin_l, ypix), (margin_l + plot_w, ypix)], fill=(220, 220, 220), width=1)
        draw.text((margin_l - 70, ypix - 8), f"{y_tick:.2f}", fill="black", font=font)
    draw.rectangle([margin_l, margin_t, margin_l + plot_w, margin_t + plot_h], outline="black", width=3)

    colors = {20: "#CD661D", 30: "#228B22", 40: "#B22222", 50: "#68228B"}
    rgb_colors = {20: (205, 102, 29), 30: (34, 139, 34), 40: (178, 34, 34), 50: (104, 34, 139)}
    for radius in r_values:
        sub = results[results["R_km"] == radius].sort_values("alpha")
        pts = [(xp(float(a)), yp(float(r))) for a, r in zip(sub["alpha"], sub["RMSE_mm_h"])]
        draw.line(pts, fill=rgb_colors[radius], width=6)
        for px, py in pts:
            draw.ellipse([px - 9, py - 9, px + 9, py + 9], fill=rgb_colors[radius], outline="black", width=2)

    optimum = results.loc[results["RMSE_mm_h"].idxmin()]
    ox, oy = xp(float(optimum["alpha"])), yp(float(optimum["RMSE_mm_h"]))
    draw.ellipse([ox - 20, oy - 20, ox + 20, oy + 20], fill="red", outline="black", width=3)
    draw.line([(ox + 18, oy - 18), (xp(0.76), yp(float(optimum["RMSE_mm_h"]) + 0.04))], fill="black", width=3)
    draw.text((xp(0.77), yp(float(optimum["RMSE_mm_h"]) + 0.05) - 20), f"Global Optimum\nR={int(optimum['R_km'])} km, alpha={optimum['alpha']:.1f}", fill="black", font=font)

    draw.text((margin_l + plot_w // 2 - 120, height - 95), "Weight Reduction Coefficient (alpha)", fill="black", font=font)
    draw.text((35, margin_t + plot_h // 2), "RMSE (mm h-1)", fill="black", font=font)
    lx, ly = margin_l + 20, margin_t + plot_h - 150
    for i, radius in enumerate(r_values):
        y0 = ly + i * 30
        draw.line([(lx, y0), (lx + 45, y0)], fill=rgb_colors[radius], width=6)
        draw.ellipse([lx + 17, y0 - 8, lx + 33, y0 + 8], fill=rgb_colors[radius], outline="black", width=1)
        draw.text((lx + 60, y0 - 8), f"R = {radius} km", fill="black", font=font)
    draw.ellipse([lx + 18, ly + 4 * 30 - 10, lx + 38, ly + 4 * 30 + 10], fill="red", outline="black", width=2)
    draw.text((lx + 60, ly + 4 * 30 - 8), "Global Optimum", fill="black", font=font)

    img.save(out_dir / "fig5_sensitivity.png")

    print(results.to_string(index=False, formatters={"RMSE_mm_h": "{:.3f}".format}))
    print(f"Saved: {out_dir / 'fig5_sensitivity.png'}")


if __name__ == "__main__":
    main()
