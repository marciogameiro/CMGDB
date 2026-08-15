# CMGDB version performance comparison

Generated 2026-08-15 20:25 on macOS-10.16-x86_64-i386-64bit.

**old** = git revision `857ec8b`; **new** = current working tree; **torch** = new with the map evaluated by torch (float64, CPU).

Single repetition per cell, flavors back-to-back per scenario for thermal fairness. `n/a` = the mechanism does not exist in that version, or (torch column) the scenario has no torch-evaluable map: the C++ interval map and the data-defined maps, which evaluate no function. In the torch column, scalar scenarios run with a torch-backed batched map attached; the Conley index phase always evaluates the scalar map. Torch runs in float64 so every flavor computes identical Morse sets; float32 evaluation is roughly twice as fast again but can perturb box covers.

| Scenario | Old | New | Torch | Old/New | New/Torch |
|---|---|---|---|---|---|
| leslie2d_python | 3.47s | 1.80s | 0.58s | 1.93x | 3.10x |
| leslie3d_python | 3.54s | 1.92s | 0.48s | 1.84x | 3.99x |
| leslie4d_python | 10.1s | 5.30s | 1.19s | 1.91x | 4.44x |
| leslie2d_python_deep | 33.5s | 17.6s | 5.78s | 1.90x | 3.04x |
| leslie3d_python_deep | 9.66s | 5.21s | 1.36s | 1.85x | 3.82x |
| leslie2d_python_batch | n/a | 0.59s | 0.59s |  | 1.00x |
| leslie3d_python_batch | n/a | 0.49s | 0.49s |  | 1.02x |
| leslie2d_python_deep_batch | n/a | 5.80s | 5.83s |  | 0.99x |
| leslie2d_interval | 2.00s | 1.34s | n/a | 1.50x |  |
| leslie2d_conley | 5.14s | 3.44s | 2.15s | 1.50x | 1.60x |
| leslie3d_conley | 15.2s | 13.9s | 13.7s | 1.09x | 1.02x |
| leslie4d_conley | 12.2s | 11.8s | 11.0s | 1.03x | 1.08x |
| leslie2d_data | 3.92s | 2.00s | n/a | 1.96x |  |
| leslie2d_data_large | 113.9s | 4.44s | n/a | 25.67x |  |
| leslie2d_data_large_batch | n/a | 4.08s | n/a |  |  |
| leslie2d_data_conley | 4.53s | 2.59s | n/a | 1.75x |  |
| product2d_small | 0.01s | 0.01s | 0.01s | 1.29x | 0.41x |
| product2d_small_conley | 0.01s | 0.01s | 0.02s | 1.20x | 0.51x |
| product2d_medium | 0.01s | 0.01s | 0.02s | 1.42x | 0.57x |
| product2d_medium_batch | n/a | 0.01s | 0.02s |  | 0.69x |
| product2d_adaptive | 0.02s | 0.01s | 0.03s | 1.48x | 0.56x |
| product2d_uniform | 0.08s | 0.04s | 0.02s | 1.88x | 2.27x |
| product2d_conley | 0.02s | 0.02s | 0.03s | 0.94x | 0.82x |
| product3d_uniform | 1.53s | 0.57s | 0.18s | 2.69x | 3.25x |
| product3d_uniform_deep | 14.3s | 4.63s | 1.38s | 3.08x | 3.37x |
| product3d_deep_adaptive | 3.3m | 73.7s | 26.2s | 2.73x | 2.82x |
| product4d_adaptive | 7.53s | 0.72s | 0.54s | 10.52x | 1.32x |
| product4d_deep | 5.2m | 31.3s | 8.85s | 10.05x | 3.53x |
| henon3d | 11.9s | 6.57s | 4.10s | 1.82x | 1.60x |
| henon3d_conley | 1.51s | 0.96s | 0.71s | 1.57x | 1.35x |
| leslie2d_expensive_live | n/a | 7.11s | 0.63s |  | 11.20x |
| leslie2d_expensive_pre | n/a | 3.18s | 0.53s |  | 6.00x |

All flavors computed the same number of Morse sets on every scenario they could run.
