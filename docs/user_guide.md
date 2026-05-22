# User Guide

This guide describes the expected inputs, outputs, command-line options, and behaviour of the reproducibility scripts.

## Input Data Format

The default input file is:

```text
data/example_rainfall_stations.csv
```

The CSV file must contain the following columns:

| Column | Type | Unit | Description |
| --- | --- | --- | --- |
| `station_id` | string | none | An anonymized station identifier. |
| `x_km` | numeric | km | Relative projected x coordinate. |
| `y_km` | numeric | km | Relative projected y coordinate. |
| `rain_mm_h` | numeric | mm h-1 | Hourly rainfall intensity. |

Coordinates should be in a projected or local Cartesian coordinate system measured in kilometers. Longitude and latitude in degrees should be projected before use because the method applies Euclidean distances.

## Core Behaviour

The default improved IDW workflow uses the following logic:

1. Compute a local z-score for each station using its `N = 5` nearest neighbors.
2. Mark a station as anomalous when the local z-score exceeds `T = 1.5`.
3. During interpolation, apply the standard IDW weight to ordinary stations.
4. For anomalous stations, set the weight to zero when `d >= R`.
5. For anomalous stations, multiply the weight by `alpha` when `d < R`.

The default parameters are `R = 30 km`, `alpha = 0.6`, and IDW power `p = 2`.

## Scripts

### `scripts/run_all.py`

Runs the complete default reproducibility workflow:

```bash
python scripts/run_all.py
```

Expected outputs:

- `outputs/table2_loocv_predictions.csv`
- `outputs/table2_loocv_summary.csv`
- `outputs/fig5_sensitivity_results.csv`
- `outputs/fig5_sensitivity.png`
- `outputs/fig6_interpolation_maps.png`
- `outputs/fig7_profile_data.csv`
- `outputs/fig7_profile.png`

### `scripts/reproduce_table2.py`

Runs leave-one-out cross-validation (LOOCV) for the interpolation methods.

```bash
python scripts/reproduce_table2.py
```

Options:

- `--data PATH`: input station CSV. Default: `data/example_rainfall_stations.csv`.
- `--out-dir PATH`: output directory. Default: `outputs`.
- `--include-geostatistical`: also run Ordinary Kriging and thin-plate spline baselines.

Expected outputs:

- `table2_loocv_predictions.csv`: station-level observed and predicted values.
- `table2_loocv_summary.csv`: RMSE, MAE, and RMSE reduction for each method.

### `scripts/reproduce_fig5_sensitivity.py`

Runs the sensitivity analysis for the spatial influence radius `R` and weight reduction coefficient `alpha`.

```bash
python scripts/reproduce_fig5_sensitivity.py
```

Options:

- `--data PATH`: input station CSV.
- `--out-dir PATH`: output directory.

Expected outputs:

- `fig5_sensitivity_results.csv`: RMSE and MAE for each `(R, alpha)` combination.
- `fig5_sensitivity.png`: sensitivity curve figure.

### `scripts/reproduce_fig6_maps.py`

Creates map-style interpolation panels for method comparison.

```bash
python scripts/reproduce_fig6_maps.py
```

Options:

- `--data PATH`: input station CSV.
- `--out-dir PATH`: output directory.
- `--grid-size INTEGER`: number of grid nodes along each axis. Default: `130`.
- `--include-kriging`: include an optional Ordinary Kriging panel.
- `--include-spline`: include an optional thin-plate spline panel.

Expected output:

- `fig6_interpolation_maps.png`

### `scripts/reproduce_fig7_profile.py`

Creates a one-dimensional profile comparison starting from the strongest extreme-value station.

```bash
python scripts/reproduce_fig7_profile.py
```

Options:

- `--data PATH`: input station CSV.
- `--out-dir PATH`: output directory.
- `--max-distance-km FLOAT`: maximum profile distance from the extreme station. Default: `80.0`.

Expected outputs:

- `fig7_profile_data.csv`
- `fig7_profile.png`

## Using Custom Data

To use a custom station dataset, prepare a CSV with the required columns. The default `run_all.py` script is intended for the bundled example dataset. For custom data, run the individual scripts and pass the file path explicitly:

```bash
python scripts/reproduce_table2.py --data path/to/your_stations.csv
python scripts/reproduce_fig5_sensitivity.py --data path/to/your_stations.csv
python scripts/reproduce_fig6_maps.py --data path/to/your_stations.csv
python scripts/reproduce_fig7_profile.py --data path/to/your_stations.csv
```

The expected behaviour is that all predictions are non-negative rainfall intensities and all generated files are written to the selected output directory. If a custom dataset has very few stations, LOOCV and adaptive-neighbor calculations may be less stable; in that case, users should adjust the neighbor-related settings in `src/dual_constraint_idw.py`.

## Notes on Optional Baselines

The default workflow focuses on the IDW-related methods because these results are deterministic and lightweight. Ordinary Kriging and thin-plate spline baselines are available through optional flags, but they may depend on solver behaviour, variogram assumptions, and platform-specific numerical details.
