"""Reproduce a one-dimensional profile comparison across an extreme station."""

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
    standard_idw_predict,
    values_from_df,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default=ROOT / "data" / "example_rainfall_stations.csv")
    parser.add_argument("--out-dir", default=ROOT / "outputs")
    parser.add_argument("--max-distance-km", type=float, default=80.0)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    config = IDWConfig(radius_km=30.0, alpha=0.6)
    df = load_station_data(args.data)
    coords = coordinates_from_df(df)
    values = values_from_df(df)
    dist_matrix = pairwise_distances(coords)
    anomaly_flags, _ = identify_anomalies(values, dist_matrix, config.anomaly_neighbors, config.anomaly_threshold)

    center_idx = int(np.argmax(values))
    center = coords[center_idx]
    # Use the direction from the domain center to the strongest extreme station.
    domain_center = np.mean(coords, axis=0)
    direction = center - domain_center
    norm = float(np.linalg.norm(direction))
    if norm == 0:
        direction = np.array([1.0, 0.0])
    else:
        direction = direction / norm

    distances = np.linspace(0.0, args.max_distance_km, 180)
    profile_points = center[None, :] + distances[:, None] * direction[None, :]

    standard_vals = []
    improved_vals = []
    for point in profile_points:
        standard_vals.append(standard_idw_predict(coords, values, point, config))
        improved_vals.append(improved_idw_predict(coords, values, anomaly_flags, point, config))

    profile = pd.DataFrame({
        "distance_km": distances,
        "standard_idw_mm_h": standard_vals,
        "improved_idw_mm_h": improved_vals,
    })
    profile.to_csv(out_dir / "fig7_profile_data.csv", index=False)

    width, height = 1500, 960
    ml, mr, mt, mb = 180, 80, 70, 150
    pw, ph = width - ml - mr, height - mt - mb
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()

    y_values = np.array(standard_vals + improved_vals, dtype=float)
    x_min, x_max = 0.0, float(args.max_distance_km)
    y_min, y_max = 0.0, float(np.nanmax(y_values) * 1.08)

    def xp(x):
        return int(ml + (x - x_min) / (x_max - x_min) * pw)

    def yp(y):
        return int(mt + (y_max - y) / (y_max - y_min) * ph)

    for x_tick in np.linspace(0, x_max, 5):
        px = xp(float(x_tick))
        draw.line([(px, mt), (px, mt + ph)], fill=(225, 225, 225), width=1)
        draw.text((px - 12, mt + ph + 18), f"{x_tick:.0f}", fill="black", font=font)
    for y_tick in np.linspace(0, y_max, 6):
        py = yp(float(y_tick))
        draw.line([(ml, py), (ml + pw, py)], fill=(225, 225, 225), width=1)
        draw.text((ml - 70, py - 8), f"{y_tick:.1f}", fill="black", font=font)
    draw.rectangle([ml, mt, ml + pw, mt + ph], outline="black", width=3)

    std_pts = [(xp(float(x)), yp(float(y))) for x, y in zip(distances, standard_vals)]
    imp_pts = [(xp(float(x)), yp(float(y))) for x, y in zip(distances, improved_vals)]
    draw.line(std_pts, fill=(198, 93, 40), width=6)
    draw.line(imp_pts, fill=(27, 127, 58), width=6)
    r_x = xp(config.radius_km)
    for yy in range(mt, mt + ph, 18):
        draw.line([(r_x, yy), (r_x, yy + 9)], fill="black", width=2)

    draw.text((ml + pw // 2 - 165, height - 92), "Distance from extreme-value station (km)", fill="black", font=font)
    draw.text((25, mt + ph // 2), "Interpolated rainfall intensity (mm h-1)", fill="black", font=font)
    lx, ly = ml + 30, mt + 30
    draw.line([(lx, ly), (lx + 55, ly)], fill=(198, 93, 40), width=6)
    draw.text((lx + 70, ly - 8), "Standard IDW", fill="black", font=font)
    draw.line([(lx, ly + 32), (lx + 55, ly + 32)], fill=(27, 127, 58), width=6)
    draw.text((lx + 70, ly + 24), "Improved IDW", fill="black", font=font)
    draw.line([(lx, ly + 64), (lx + 55, ly + 64)], fill="black", width=2)
    draw.text((lx + 70, ly + 56), "R = 30 km", fill="black", font=font)

    img.save(out_dir / "fig7_profile.png")
    print(f"Saved: {out_dir / 'fig7_profile.png'}")


if __name__ == "__main__":
    main()
