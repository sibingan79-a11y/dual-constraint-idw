"""Core routines for dual-constraint IDW interpolation.

The method follows the manuscript logic:
1. identify anomalous extreme-value stations using a local z-score test;
2. set anomalous-station weights to zero when d >= R;
3. multiply anomalous-station weights by alpha when d < R.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist


@dataclass(frozen=True)
class IDWConfig:
    power: float = 2.0
    search_points: int = 12
    anomaly_neighbors: int = 5
    anomaly_threshold: float = 1.5
    radius_km: float = 30.0
    alpha: float = 0.6
    local_radius_km: float = 30.0
    local_min_points: int = 3
    aidw_neighbors: int = 5
    aidw_p_min: float = 1.0
    aidw_p_max: float = 4.0


def load_station_data(path: str | Path) -> pd.DataFrame:
    """Load the processed station CSV used by the reproducibility scripts."""
    df = pd.read_csv(path)
    required = {"station_id", "x_km", "y_km", "rain_mm_h"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    return df.copy()


def coordinates_from_df(df: pd.DataFrame) -> np.ndarray:
    return df[["x_km", "y_km"]].to_numpy(dtype=float)


def values_from_df(df: pd.DataFrame) -> np.ndarray:
    return df["rain_mm_h"].to_numpy(dtype=float)


def pairwise_distances(coords: np.ndarray) -> np.ndarray:
    return cdist(coords, coords, metric="euclidean")


def identify_anomalies(
    values: np.ndarray,
    dist_matrix: np.ndarray,
    n_neighbors: int = 5,
    threshold: float = 1.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Identify local anomalous extreme-value stations.

    A station is anomalous when its local z-score relative to its N nearest
    neighbors exceeds the threshold T.
    """
    n = len(values)
    flags = np.zeros(n, dtype=bool)
    scores = np.zeros(n, dtype=float)

    for i in range(n):
        order = np.argsort(dist_matrix[i])
        neighbors = order[order != i][:n_neighbors]
        neighbor_values = values[neighbors]
        std = float(np.std(neighbor_values))
        score = (values[i] - float(np.mean(neighbor_values))) / std if std > 0 else 0.0
        scores[i] = score
        flags[i] = score > threshold

    return flags, scores


def idw_weights(distances: np.ndarray, power: float) -> np.ndarray:
    safe_dist = np.maximum(distances.astype(float), 1e-12)
    return 1.0 / (safe_dist**power)


def weighted_average(values: np.ndarray, weights: np.ndarray) -> float:
    total = float(np.sum(weights))
    if total <= 0 or not np.isfinite(total):
        return 0.0
    return float(np.sum(weights * values) / total)


def nearest_subset(
    train_coords: np.ndarray,
    target_coord: np.ndarray,
    search_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    tree = cKDTree(train_coords)
    k = min(search_points, len(train_coords))
    distances, indices = tree.query(target_coord, k=k)
    return np.atleast_1d(indices), np.atleast_1d(distances)


def standard_idw_predict(
    train_coords: np.ndarray,
    train_values: np.ndarray,
    target_coord: np.ndarray,
    config: IDWConfig,
) -> float:
    idx, distances = nearest_subset(train_coords, target_coord, config.search_points)
    weights = idw_weights(distances, config.power)
    return weighted_average(train_values[idx], weights)


def improved_idw_predict(
    train_coords: np.ndarray,
    train_values: np.ndarray,
    train_anomaly_flags: np.ndarray,
    target_coord: np.ndarray,
    config: IDWConfig,
) -> float:
    idx, distances = nearest_subset(train_coords, target_coord, config.search_points)
    weights = idw_weights(distances, config.power)
    anomalous = train_anomaly_flags[idx]

    weights[anomalous & (distances >= config.radius_km)] = 0.0
    weights[anomalous & (distances < config.radius_km)] *= config.alpha
    return weighted_average(train_values[idx], weights)


def local_radius_idw_predict(
    train_coords: np.ndarray,
    train_values: np.ndarray,
    target_coord: np.ndarray,
    config: IDWConfig,
) -> float:
    distances = np.linalg.norm(train_coords - target_coord, axis=1)
    inside = distances <= config.local_radius_km
    if int(np.sum(inside)) >= config.local_min_points:
        use_values = train_values[inside]
        use_distances = distances[inside]
    else:
        order = np.argsort(distances)[: min(config.search_points, len(distances))]
        use_values = train_values[order]
        use_distances = distances[order]
    return weighted_average(use_values, idw_weights(use_distances, config.power))


def aidw_power(
    train_dist_matrix: np.ndarray,
    test_distances: np.ndarray,
    config: IDWConfig,
) -> float:
    """A simple adaptive-power IDW baseline based on nearest-neighbor spacing."""
    train_n = train_dist_matrix.shape[0]
    k = min(config.aidw_neighbors, train_n - 1)
    if k <= 0:
        return config.power

    local_mean = float(np.mean(np.sort(test_distances)[:k]))
    masked = train_dist_matrix + np.eye(train_n) * 1e12
    global_mean = float(np.mean(np.sort(masked, axis=1)[:, :k]))
    if global_mean <= 0 or not np.isfinite(global_mean):
        return config.power

    density_ratio = local_mean / global_mean
    return float(np.clip(config.power * density_ratio, config.aidw_p_min, config.aidw_p_max))


def aidw_predict(
    train_coords: np.ndarray,
    train_values: np.ndarray,
    target_coord: np.ndarray,
    config: IDWConfig,
) -> tuple[float, float]:
    distances = np.linalg.norm(train_coords - target_coord, axis=1)
    order = np.argsort(distances)[: min(config.search_points, len(distances))]
    train_dist_matrix = pairwise_distances(train_coords)
    p_adaptive = aidw_power(train_dist_matrix, distances, config)
    weights = idw_weights(distances[order], p_adaptive)
    return weighted_average(train_values[order], weights), p_adaptive


def ordinary_kriging_predict(
    train_coords: np.ndarray,
    train_values: np.ndarray,
    target_coord: np.ndarray,
) -> float:
    from pykrige.ok import OrdinaryKriging

    ok = OrdinaryKriging(
        train_coords[:, 0],
        train_coords[:, 1],
        train_values,
        variogram_model="linear",
        verbose=False,
        enable_plotting=False,
    )
    pred, _ = ok.execute("points", np.array([target_coord[0]]), np.array([target_coord[1]]))
    return max(0.0, float(pred[0]))


def spline_predict(
    train_coords: np.ndarray,
    train_values: np.ndarray,
    target_coord: np.ndarray,
) -> float:
    rbf = Rbf(train_coords[:, 0], train_coords[:, 1], train_values, function="thin_plate")
    return max(0.0, float(rbf(target_coord[0], target_coord[1])))


def loocv_predictions(
    coords: np.ndarray,
    values: np.ndarray,
    anomaly_flags: np.ndarray,
    config: IDWConfig,
    include_geostatistical: bool = True,
) -> pd.DataFrame:
    """Run LOOCV for the standard, improved, and baseline interpolation methods."""
    rows: list[dict[str, float | int | bool]] = []
    n = len(values)

    for i in range(n):
        train_idx = np.delete(np.arange(n), i)
        train_coords = coords[train_idx]
        train_values = values[train_idx]
        train_flags = anomaly_flags[train_idx]
        target = coords[i]

        standard = standard_idw_predict(train_coords, train_values, target, config)
        improved = improved_idw_predict(train_coords, train_values, train_flags, target, config)
        local_radius = local_radius_idw_predict(train_coords, train_values, target, config)
        aidw, p_adaptive = aidw_predict(train_coords, train_values, target, config)

        row: dict[str, float | int | bool] = {
            "station_index": i,
            "observed": float(values[i]),
            "is_anomaly": bool(anomaly_flags[i]),
            "Standard IDW": standard,
            "Improved IDW": improved,
            "Local-radius IDW": local_radius,
            "AIDW": aidw,
            "AIDW_power": p_adaptive,
        }

        if include_geostatistical:
            try:
                row["Ordinary Kriging"] = ordinary_kriging_predict(train_coords, train_values, target)
            except Exception:
                row["Ordinary Kriging"] = np.nan
            try:
                row["Spline (Thin Plate)"] = spline_predict(train_coords, train_values, target)
            except Exception:
                row["Spline (Thin Plate)"] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def summarize_errors(df: pd.DataFrame, methods: Iterable[str]) -> pd.DataFrame:
    rows = []
    baseline_rmse = None
    for method in methods:
        valid = df.dropna(subset=[method])
        errors = valid[method].to_numpy(dtype=float) - valid["observed"].to_numpy(dtype=float)
        rmse = float(np.sqrt(np.mean(errors**2)))
        mae = float(np.mean(np.abs(errors)))
        if baseline_rmse is None:
            baseline_rmse = rmse
        reduction = (baseline_rmse - rmse) / baseline_rmse * 100.0 if baseline_rmse else 0.0
        rows.append(
            {
                "Method": method,
                "RMSE_mm_h": rmse,
                "MAE_mm_h": mae,
                "RMSE_reduction_percent": reduction,
            }
        )
    return pd.DataFrame(rows)

