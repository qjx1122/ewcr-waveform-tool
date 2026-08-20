# Reference result provenance

The `reference_*.csv` files and `reference_*.png` figures are included so that the reorganized project can be inspected before a full MATLAB rerun.

- ASBC and DWT-Hybrid fixed-scenario values, the original ASBC/DWT/CS rate-distortion points, and the real-data ASBC/DWT/CS baselines are read from the user-provided `results/benchmark_results.mat`; an unchanged copy is stored in `source/benchmark_results_original.mat`.
- The full 10-cycle CS-OMP fixed-scenario rows were numerically mirrored from the supplied MATLAB equations/configuration to remove the old 4-cycle truncation from the comparison table.
- MMC and SVDCS reference rows were numerically mirrored from the supplied MATLAB implementations for pre-validation because MATLAB/Octave is not installed in the artifact-generation runtime.
- The authoritative final five-codec comparison should be regenerated in MATLAB by running `run_all_experiments.m`; this writes fresh CSV/MAT files under `results/` and figures under `figures/`.

The PDF report labels this distinction explicitly and does not present the mirrored rows as a bit-identical MATLAB execution.
