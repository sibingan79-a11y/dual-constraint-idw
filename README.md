# Dual-Constraint IDW for Mitigating the Bull's-Eye Effect

This repository contains the reproducibility code for the manuscript:

**A Dual-Constraint Method Based on Local Influence Distance and Weight Dominance for Mitigating the Bull's-Eye Effect in IDW Interpolation**

The proposed method improves inverse distance weighting (IDW) interpolation for strongly non-stationary spatial fields containing isolated or semi-isolated extreme observations. It first identifies anomalous extreme-value stations and then applies two asymmetric constraints:

- **Far-field truncation:** anomalous-station weights are set to zero when `d >= R`.
- **Near-field attenuation:** anomalous-station weights are multiplied by `alpha` when `d < R`.

## Repository Structure

```text
dual_constraint_idw_code/
├─ README.md
├─ LICENSE
├─ requirements.txt
├─ environment.yml
├─ data/
│  ├─ example_rainfall_stations.csv
│  └─ README.md
├─ src/
│  ├─ __init__.py
│  └─ dual_constraint_idw.py
├─ scripts/
│  ├─ run_all.py
│  ├─ reproduce_table2.py
│  ├─ reproduce_fig5_sensitivity.py
│  ├─ reproduce_fig6_maps.py
│  └─ reproduce_fig7_profile.py
├─ outputs/
│  └─ README.md
└─ docs/
   └─ method_notes.md
```

## Environment

Python 3.10 is recommended.

Install dependencies with pip:

```bash
pip install -r requirements.txt
```

Alternatively, create a conda environment:

```bash
conda env create -f environment.yml
conda activate dual-constraint-idw
```

## Reproduce the Main Results

Run all scripts from the repository root:

```bash
python scripts/run_all.py
```

Or run individual scripts:

```bash
python scripts/reproduce_table2.py
python scripts/reproduce_fig5_sensitivity.py
python scripts/reproduce_fig6_maps.py
python scripts/reproduce_fig7_profile.py
```

The generated files will be written to `outputs/`.

By default, `run_all.py` focuses on the stable and fully reproducible IDW-related results. Optional Ordinary Kriging and thin-plate spline baselines can be enabled in the relevant scripts with command-line flags, but they depend on external numerical solvers and may behave differently across platforms.

## Input Data

The example dataset `data/example_rainfall_stations.csv` is a processed and anonymized station dataset. Station names and original identifiers are removed; coordinates are converted to relative kilometer coordinates to demonstrate the algorithm and reproduce the workflow.

The original meteorological observations were obtained from the National Meteorological Information Center of China and may be subject to the data provider's access policy.

## Main Parameters

The default manuscript parameters are:

- local anomaly detection: `N = 5`, `T = 1.5`;
- spatial influence radius: `R = 30 km`;
- weight reduction coefficient: `alpha = 0.6`;
- IDW power exponent: `p = 2`;
- nearest neighbors for interpolation: `12`.

## Reproduced Outputs

The default one-command workflow reproduces:

- LOOCV metrics for Standard IDW, Improved IDW, Local-radius IDW, and AIDW;
- the R-alpha sensitivity analysis corresponding to the parameter-coupling figure;
- spatial interpolation maps for IDW variants;
- a one-dimensional profile comparison between Standard IDW and Improved IDW.

## License

The code is released under the MIT License.
