# Method Notes

The dual-constraint IDW method modifies only the contribution of locally anomalous extreme-value stations.

For a target location and a training station at distance `d`, the improved weight is:

```text
w_i* = 0                 if station i is anomalous and d >= R
w_i* = alpha * d_i^(-p)  if station i is anomalous and d < R
w_i* = d_i^(-p)          otherwise
```

The anomaly flag is computed from the local z-score:

```text
E_i = (z_i - mean(z_neighbors)) / std(z_neighbors)
```

A station is treated as anomalous when `E_i > T`. The default values used in the manuscript are `N = 5` and `T = 1.5`.

