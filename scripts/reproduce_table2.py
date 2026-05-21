"""Reproduce the LOOCV comparison table."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.dual_constraint_idw import (  # noqa: E402
    IDWConfig,
    coordinates_from_df,
    identify_anomalies,
    load_station_data,
    loocv_predictions,
    pairwise_distances,
    summarize_errors,
    values_from_df,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default=ROOT / "data" / "example_rainfall_stations.csv")
    parser.add_argument("--out-dir", default=ROOT / "outputs")
    parser.add_argument(
        "--include-geostatistical",
        action="store_true",
        help="Also run Ordinary Kriging and thin-plate spline baselines. "
        "These optional baselines depend on external numerical solvers.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    config = IDWConfig(radius_km=30.0, alpha=0.6)
    df = load_station_data(args.data)
    coords = coordinates_from_df(df)
    values = values_from_df(df)
    dist_matrix = pairwise_distances(coords)
    anomaly_flags, anomaly_scores = identify_anomalies(
        values,
        dist_matrix,
        n_neighbors=config.anomaly_neighbors,
        threshold=config.anomaly_threshold,
    )

    predictions = loocv_predictions(
        coords,
        values,
        anomaly_flags,
        config,
        include_geostatistical=args.include_geostatistical,
    )
    predictions["anomaly_score"] = anomaly_scores
    predictions.to_csv(out_dir / "table2_loocv_predictions.csv", index=False)

    methods = [
        "Standard IDW",
        "Improved IDW",
        "Local-radius IDW",
        "AIDW",
    ]
    if args.include_geostatistical:
        methods[2:2] = ["Ordinary Kriging", "Spline (Thin Plate)"]
    summary = summarize_errors(predictions, methods)
    summary.to_csv(out_dir / "table2_loocv_summary.csv", index=False)

    print(summary.to_string(index=False, formatters={
        "RMSE_mm_h": "{:.3f}".format,
        "MAE_mm_h": "{:.3f}".format,
        "RMSE_reduction_percent": "{:.1f}".format,
    }))
    print(f"Saved: {out_dir / 'table2_loocv_summary.csv'}")


if __name__ == "__main__":
    main()
