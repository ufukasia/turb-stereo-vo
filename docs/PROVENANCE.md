# Provenance

- Base: PyEco-SLAM (https://github.com/ufukasia/PyEco), commit
  `52fe0173c053ac8d9af50dbe4bdab871a3de0ca4` (2026-11-04). Copied without results, paper, docs,
  outputs, and git history. The upstream repository is not modified by this work.
- `turb/`: turbulence simulator, covariance arms, estimators, metrics, and runners written for this
  study (2026-11-07 onward; includes the point-spread-function scale correction,
  `psf_legacy_scale=false`).
- Rule kept throughout: with turbulence switched off, the output of this tree is bit-identical to
  PyEco-SLAM (gate E0, `turb/e0_bit_identity.py`).
- Reason for the platform: on a rendered oblique AirSim sequence (81 m) the SE(3) ATE of the earlier
  pipeline was 2.49 m and that of PyEco-SLAM 0.035 m.
- Vendored third-party code and licenses: see `NOTICE`.
- `airsim/`: recording and grid scripts for the rendered test bed (no binaries, no maps).
- `results/`: records and analysis outputs of every experiment; `evidence` is a symbolic link to
  `results` because the analysis scripts read `evidence/<EXPERIMENT>/`.
