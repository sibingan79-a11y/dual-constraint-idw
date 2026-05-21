"""Reproduce map-style spatial interpolation panels for method comparison."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.dual_constraint_idw import (  # noqa: E402
    IDWConfig,
    coordinates_from_df,
    identify_anomalies,
    idw_weights,
    load_station_data,
    pairwise_distances,
    values_from_df,
    weighted_average,
)


def grid_idw(coords, values, anomaly_flags, grid_points, config, improved=False):
    from scipy.spatial import cKDTree

    tree = cKDTree(coords)
    k = min(config.search_points, len(values))
    distances, indices = tree.query(grid_points, k=k)
    distances = np.atleast_2d(distances)
    indices = np.atleast_2d(indices)
    if distances.shape[0] != len(grid_points):
        distances = distances.T
        indices = indices.T

    weights = idw_weights(distances, config.power)
    if improved:
        anomalous = anomaly_flags[indices]
        weights[anomalous & (distances >= config.radius_km)] = 0.0
        weights[anomalous & (distances < config.radius_km)] *= config.alpha

    sums = np.sum(weights, axis=1)
    sums[sums <= 0] = np.nan
    pred = np.sum(weights * values[indices], axis=1) / sums
    return np.nan_to_num(pred, nan=0.0)


def grid_local_radius_idw(coords, values, grid_points, config):
    output = []
    for point in grid_points:
        distances = np.linalg.norm(coords - point, axis=1)
        inside = distances <= config.local_radius_km
        if int(np.sum(inside)) >= config.local_min_points:
            use_values = values[inside]
            use_distances = distances[inside]
        else:
            order = np.argsort(distances)[: min(config.search_points, len(distances))]
            use_values = values[order]
            use_distances = distances[order]
        output.append(weighted_average(use_values, idw_weights(use_distances, config.power)))
    return np.asarray(output, dtype=float)


def grid_aidw(coords, values, grid_points, config):
    train_dist = pairwise_distances(coords)
    sorted_train = np.sort(train_dist + np.eye(len(coords)) * 1e12, axis=1)
    k = min(config.aidw_neighbors, len(coords) - 1)
    global_mean = float(np.mean(sorted_train[:, :k]))
    output = []
    for point in grid_points:
        distances = np.linalg.norm(coords - point, axis=1)
        order = np.argsort(distances)[: min(config.search_points, len(distances))]
        local_mean = float(np.mean(np.sort(distances)[:k]))
        ratio = local_mean / global_mean if global_mean > 0 else 1.0
        power = float(np.clip(config.power * ratio, config.aidw_p_min, config.aidw_p_max))
        output.append(weighted_average(values[order], idw_weights(distances[order], power)))
    return np.asarray(output, dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default=ROOT / "data" / "example_rainfall_stations.csv")
    parser.add_argument("--out-dir", default=ROOT / "outputs")
    parser.add_argument("--grid-size", type=int, default=130)
    parser.add_argument(
        "--include-kriging",
        action="store_true",
        help="Run the optional Ordinary Kriging panel. Disabled by default for robustness.",
    )
    parser.add_argument(
        "--include-spline",
        action="store_true",
        help="Run the optional thin-plate spline panel. Disabled by default for robustness.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    config = IDWConfig(radius_km=30.0, alpha=0.6)
    df = load_station_data(args.data)
    coords = coordinates_from_df(df)
    values = values_from_df(df)
    dist_matrix = pairwise_distances(coords)
    anomaly_flags, _ = identify_anomalies(values, dist_matrix, config.anomaly_neighbors, config.anomaly_threshold)

    pad = 10.0
    x = np.linspace(coords[:, 0].min() - pad, coords[:, 0].max() + pad, args.grid_size)
    y = np.linspace(coords[:, 1].min() - pad, coords[:, 1].max() + pad, args.grid_size)
    xx, yy = np.meshgrid(x, y)
    grid_points = np.column_stack([xx.ravel(), yy.ravel()])

    standard = grid_idw(coords, values, anomaly_flags, grid_points, config, improved=False).reshape(xx.shape)
    improved = grid_idw(coords, values, anomaly_flags, grid_points, config, improved=True).reshape(xx.shape)
    local_radius = grid_local_radius_idw(coords, values, grid_points, config).reshape(xx.shape)
    aidw = grid_aidw(coords, values, grid_points, config).reshape(xx.shape)

    if args.include_kriging:
        from pykrige.ok import OrdinaryKriging

        ok = OrdinaryKriging(coords[:, 0], coords[:, 1], values, variogram_model="linear")
        kriging, _ = ok.execute("grid", x, y)
        kriging = np.maximum(np.asarray(kriging, dtype=float), 0.0)
    else:
        kriging = np.full_like(standard, np.nan)

    if args.include_spline:
        from scipy.interpolate import Rbf

        rbf = Rbf(coords[:, 0], coords[:, 1], values, function="thin_plate")
        spline = np.maximum(rbf(xx, yy), 0.0)
    else:
        spline = np.full_like(standard, np.nan)

    panels = [
        ("(a) Standard IDW", standard),
        ("(b) Local-radius IDW", local_radius),
        ("(c) AIDW", aidw),
        ("(d) Improved IDW", improved),
    ]

    vmax = np.nanpercentile(np.stack([p[1] for p in panels if np.isfinite(p[1]).any()]), 98)

    def colorize(arr):
        arr = np.nan_to_num(arr, nan=0.0)
        t = np.clip(arr / max(vmax, 1e-9), 0, 1)
        # Simple yellow-green-blue ramp.
        r = (255 * (1 - t) + 20 * t).astype(np.uint8)
        g = (245 * (1 - t) + 110 * t).astype(np.uint8)
        b = (170 * (1 - t) + 185 * t).astype(np.uint8)
        return np.dstack([r, g, b])

    panel_w, panel_h = 620, 520
    title_h = 42
    pad = 45
    img = Image.new("RGB", (panel_w * 2 + pad * 3, (panel_h + title_h) * 2 + pad * 3), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()

    x_min, x_max = float(x.min()), float(x.max())
    y_min, y_max = float(y.min()), float(y.max())

    def map_x(xval, ox):
        return int(ox + (xval - x_min) / (x_max - x_min) * panel_w)

    def map_y(yval, oy):
        return int(oy + panel_h - (yval - y_min) / (y_max - y_min) * panel_h)

    for n, (title, arr) in enumerate(panels):
        col = n % 2
        row = n // 2
        ox = pad + col * (panel_w + pad)
        oy = pad + row * (panel_h + title_h + pad) + title_h
        panel = Image.fromarray(colorize(arr)).resize((panel_w, panel_h), resample=Image.BILINEAR)
        img.paste(panel, (ox, oy))
        draw.rectangle([ox, oy, ox + panel_w, oy + panel_h], outline="black", width=2)
        draw.text((ox, oy - 30), title, fill="black", font=font)
        for cx, cy in coords:
            px, py = map_x(cx, ox), map_y(cy, oy)
            draw.ellipse([px - 4, py - 4, px + 4, py + 4], fill="white", outline="black")
        for cx, cy in coords[anomaly_flags]:
            px, py = map_x(cx, ox), map_y(cy, oy)
            draw.polygon([(px, py - 9), (px + 4, py - 2), (px + 10, py - 2), (px + 5, py + 3),
                          (px + 7, py + 10), (px, py + 6), (px - 7, py + 10), (px - 5, py + 3),
                          (px - 10, py - 2), (px - 4, py - 2)], fill="red", outline="black")

    # Color bar.
    cb_x = panel_w * 2 + pad * 2 + 12
    cb_y = pad + title_h
    cb_h = panel_h * 2 + title_h + pad
    for j in range(cb_h):
        val = 1 - j / cb_h
        rr = int(255 * (1 - val) + 20 * val)
        gg = int(245 * (1 - val) + 110 * val)
        bb = int(170 * (1 - val) + 185 * val)
        draw.line([(cb_x, cb_y + j), (cb_x + 22, cb_y + j)], fill=(rr, gg, bb))
    draw.rectangle([cb_x, cb_y, cb_x + 22, cb_y + cb_h], outline="black")
    draw.text((cb_x + 30, cb_y), f"{vmax:.1f}", fill="black", font=font)
    draw.text((cb_x + 30, cb_y + cb_h - 10), "0", fill="black", font=font)
    draw.text((cb_x - 20, cb_y + cb_h + 16), "mm h-1", fill="black", font=font)

    img.save(out_dir / "fig6_interpolation_maps.png")
    print(f"Saved: {out_dir / 'fig6_interpolation_maps.png'}")


if __name__ == "__main__":
    main()
