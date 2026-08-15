# CMGDB version performance comparison

Generated 2026-08-15 18:44 on macOS-10.16-x86_64-i386-64bit.

**old** = git revision `857ec8b`; **new** = current working tree.

Single repetition per cell, old and new back-to-back per scenario for thermal fairness. `n/a` = the mechanism does not exist in that version.

| Scenario | Old | New | Speedup |
|---|---|---|---|
| leslie2d_python | 3.61s | 1.87s | 1.93x |
| leslie3d_python | 3.46s | 1.81s | 1.91x |
| leslie4d_python | 10.1s | 5.18s | 1.95x |
| leslie2d_python_deep | 32.9s | 17.7s | 1.86x |
| leslie3d_python_deep | 9.58s | 5.00s | 1.92x |
| leslie2d_python_batch | n/a | 0.58s |  |
| leslie3d_python_batch | n/a | 0.47s |  |
| leslie2d_python_deep_batch | n/a | 5.79s |  |
| leslie2d_interval | 2.02s | 1.37s | 1.48x |
| leslie2d_conley | 5.15s | 3.48s | 1.48x |
| leslie3d_conley | 14.4s | 13.8s | 1.05x |
| leslie4d_conley | 11.9s | 11.5s | 1.04x |
| leslie2d_data | 3.60s | 1.83s | 1.97x |
| leslie2d_data_large | 2.1m | 6.39s | 19.61x |
| leslie2d_data_large_batch | n/a | 4.13s |  |
| leslie2d_data_conley | 4.23s | 2.42s | 1.75x |
| product2d_small | 0.01s | 0.01s | 1.31x |
| product2d_small_conley | 0.01s | 0.01s | 1.13x |
| product2d_medium | 0.01s | 0.01s | 1.41x |
| product2d_medium_batch | n/a | 0.01s |  |
| product2d_adaptive | 0.02s | 0.02s | 1.32x |
| product2d_uniform | 0.08s | 0.05s | 1.79x |
| product2d_conley | 0.02s | 0.02s | 1.28x |
| product3d_uniform | 1.55s | 0.56s | 2.79x |
| product3d_uniform_deep | 12.6s | 4.52s | 2.80x |
| product3d_deep_adaptive | 3.5m | 74.2s | 2.81x |
| product4d_adaptive | 7.40s | 0.72s | 10.22x |
| product4d_deep | 5.3m | 31.7s | 10.09x |
| henon3d | 11.9s | 6.58s | 1.81x |
| henon3d_conley | 1.50s | 0.99s | 1.51x |
| leslie2d_expensive_live | n/a | 1.39s |  |
| leslie2d_expensive_pre | n/a | 1.24s |  |

Both versions computed the same number of Morse sets on every scenario they could run.
