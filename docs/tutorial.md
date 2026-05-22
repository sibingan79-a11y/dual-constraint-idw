# Tutorial

This tutorial shows how to install the repository, run the example workflow, inspect the outputs, and test the method with a custom dataset.

## 1. Create the Environment

Clone or download the repository, then open a terminal in the repository root.

Using pip:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

Using conda:

```bash
conda env create -f environment.yml
conda activate dual-constraint-idw
```

## 2. Run the Default Workflow

Run:

```bash
python scripts/run_all.py
```

The script runs the LOOCV comparison, the parameter sensitivity analysis, the map-style interpolation comparison, and the one-dimensional profile comparison.

## 3. Check the LOOCV Table

Open:

```text
outputs/table2_loocv_summary.csv
```

The default example should produce the following pattern:

- Improved IDW has a lower RMSE and MAE than Standard IDW.
- Local-radius IDW may be close to Standard IDW when the fixed radius does not change the effective neighborhood structure.
- AIDW may not improve the result when adaptive distance-decay alone does not explicitly constrain anomalous extreme stations.

The exact numeric values are printed in the terminal when `scripts/reproduce_table2.py` runs.

## 4. Inspect the Parameter Sensitivity Result

Open:

```text
outputs/fig5_sensitivity_results.csv
outputs/fig5_sensitivity.png
```

The sensitivity table lists the RMSE and MAE for each tested combination of `R` and `alpha`. The default example is designed to demonstrate the parameter-coupling behaviour of the dual-constraint method.

## 5. Inspect the Spatial and Profile Outputs

Open:

```text
outputs/fig6_interpolation_maps.png
outputs/fig7_profile.png
```

The map-style panels compare spatial interpolation behaviour among IDW variants. The profile figure compares how Standard IDW and Improved IDW behave with increasing distance from the strongest extreme-value station.

## 6. Run One Script at a Time

Each analysis can also be run independently:

```bash
python scripts/reproduce_table2.py
python scripts/reproduce_fig5_sensitivity.py
python scripts/reproduce_fig6_maps.py
python scripts/reproduce_fig7_profile.py
```

For optional geostatistical baselines:

```bash
python scripts/reproduce_table2.py --include-geostatistical
```

For optional map panels:

```bash
python scripts/reproduce_fig6_maps.py --include-kriging --include-spline
```

## 7. Use a Custom Dataset

Prepare a CSV file with these columns:

```text
station_id,x_km,y_km,rain_mm_h
```

Then run:

```bash
python scripts/reproduce_table2.py --data path/to/custom_stations.csv --out-dir outputs_custom
```

The same `--data` and `--out-dir` options can be used with the figure scripts.

## 8. Interpret Common Issues

If the scripts fail because a package is missing, reinstall the environment from `requirements.txt` or `environment.yml`.

If a custom dataset uses longitude and latitude, project the coordinates to a metric coordinate system first and convert distances to kilometers.

If optional Ordinary Kriging or spline runs fail, first verify that the default IDW workflow runs successfully. The optional baselines are included for comparison and may be more sensitive to numerical solver behaviour.
