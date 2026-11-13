# turb-stereo-vo: Stereo Visual Odometry Under Atmospheric Optical Turbulence

> Code, run records, and pre-registration of the study "Stereo Visual Odometry Under Atmospheric Optical
> Turbulence: Online Image-Motion Estimation and Adaptive Outlier Gating" (Ufuk Asil).

Long-range cameras look through air whose refractive index fluctuates; the image of a static scene then
moves by several pixels, differently in every part of the field. This repository measures what that does to
stereo visual odometry (VO) and tests a remedy inside a factor-graph stereo pipeline (PyEco-SLAM):

- **Simulator.** Anisoplanatic tilt and short-exposure blur from Kolmogorov phase screens, injected after
  rectification, with a quasi-static and a fast temporal regime (`turb/`).
- **Online estimator.** The tilt standard deviation is estimated during flight from the differential motion
  of tracked features, by a structure-function fit, without a reference target (`turb/dimm_v2.py`).
- **Adaptive gate.** The estimate scales the pixel noise and the outlier gates of tracking together
  (`--turb_gate 3`); scaling the noise alone does not change the accuracy.

<p align="center"><img src="docs/figures/hero.png" width="100%"></p>
<p align="center"><sub><b>Figure 1.</b> EuRoC MH05, frames 502–901, fast regime, D/r0 = 8, seed 3: the window with the largest gain in the confirmatory grid, shown as a selected example. Top: the same left-camera frame without and with turbulence. Bottom: causal trajectories and position error of the baseline and of the adaptive gate. The median gain over the grid is smaller; see Highlights.</sub></p>

| Window (frames) | Baseline | Scalar Variant | **Adaptive Gate** | Oracle Gate | Gate / Baseline |
|:---|---:|---:|---:|---:|---:|
| 102–501 | 0.176 | 0.165 | **0.128** | 0.108 | 0.728 |
| 502–901 | 1.433 | 1.218 | **0.122** | 0.083 | 0.085 |
| 902–1301 | 1.189 | 0.912 | **0.543** | 0.373 | 0.457 |

<p align="center"><sub>ATE [m] of the three windows of the cell of Figure 1. The scalar variant uses the same online estimate without gate scaling; the oracle gate is given the true tilt strength. Cells are examples; the aggregate follows.</sub></p>

## Highlights

Median over 36 units (9 sequences × 4 seeds) of the paired ATE ratio, with 95% bootstrap confidence
intervals; values below 1 favour the first arm. The gated generic arm applies the same gate rule with a scale taken
from the running RMS of the reprojection residuals instead of the turbulence estimate. Both grids were
pre-registered and analysed once.

| Condition | **Adaptive Gate** / Baseline | Gated Generic / Baseline | Adaptive Gate / Gated Generic | Scalar Variant / Baseline | Oracle Gate / Baseline | Position ANEES, Gate / Baseline |
|:---|---:|---:|---:|---:|---:|---:|
| Quasi-static, D/r0 = 8 | **0.549 [0.415, 0.638]** | 0.542 | 0.890 [0.842, 1.001] | 1.026 | 0.412 | 0.156 |
| Quasi-static, D/r0 = 16 | **0.794 [0.752, 0.849]** | 0.849 | 0.953 [0.903, 0.988] | 1.010 | 0.582 | 0.239 |
| Fast, D/r0 = 8 | **0.362 [0.347, 0.529]** | 0.440 | 0.895 [0.824, 0.951] | 1.006 | 0.343 | 0.075 |
| Fast, D/r0 = 4 | **0.613 [0.564, 0.712]** | 0.672 | 0.952 [0.923, 1.023] | 0.937 | 0.662 | 0.401 |
| No turbulence | 1.020 [0.976, 1.047] | 0.989 | 1.019 [0.999, 1.037] | 0.991 | – | 0.683 |

- Turbulence degrades every tested system: at D/r0 = 8 the ATE grows by a factor of 2.6 to 7.0, and two of the
  four systems lose more than 40% of the windows at D/r0 = 16.
- Geometric displacement (tilt), not blur, drives the loss.
- Scaling only the measurement noise leaves the ATE unchanged (ratios 0.94 to 1.03); scaling the gates with
  the same estimate lowers it by 21% to 64%.
- Most of the gain comes from scaling the gates: the gated generic arm also lowers the ATE (ratios 0.44 to
  0.85). Compared directly, the turbulence estimate is more accurate by 5% to 11% in three conditions and
  equivalent in the fourth, and it gives a lower position ANEES in all four.
- Without turbulence the gate is neutral (upper end of the interval 1.047, pre-registered bound 1.10).
- On real imagery with scintillometer references the estimate follows the measured turbulence strength
  (rank correlation 0.863).

Limitations: the accuracy results are obtained with simulated turbulence on EuRoC and a rendered scene, the
gate does not reach the oracle at D/r0 = 16, and ten pre-registered predictions failed (reported in the
paper and in `results/`).

## Results

### 1. Degradation of four odometry systems

<p align="center"><img src="docs/figures/degradation.png" width="100%"></p>
<p align="center"><sub><b>Figure 2.</b> (a) Growth of the ATE relative to the turbulence-free run and (b) share of failed windows, against turbulence strength. The studied baseline is PyEco-SLAM; the other systems run on the same degraded images.</sub></p>

| System | D/r0 = 4 | D/r0 = 8 | D/r0 = 16 | D/r0 = 24 | Failed Windows at 16 [%] | Failed Windows at 24 [%] |
|:---|---:|---:|---:|---:|---:|---:|
| PyEco-SLAM (studied baseline) | 2.25 | 6.99 | 21.37 | 34.68 | 3.7 | 19.3 |
| ORB-SLAM3 (stereo) | 2.23 | 5.35 | 20.18 | 45.50 | 44.4 | 79.0 |
| AirSLAM | 1.66 | 2.62 | 10.49 | no valid window | 51.9 | 100.0 |
| DROID-SLAM | 1.89 | 3.92 | 7.37 | 10.93 | 0.0 | 0.0 |

<p align="center"><sub>ATE growth factor. Studied baseline: ratio of level medians over 135 windows; other systems: median paired ratio over the windows in which both runs succeeded (81 windows per level).</sub></p>

### 2. Adaptive gate

<p align="center"><img src="docs/figures/gate.png" width="100%"></p>
<p align="center"><sub><b>Figure 3.</b> (a) ATE ratio to the baseline with 95% confidence intervals and (b) position ANEES ratio to the baseline, for the scalar variant, the adaptive gate, and the oracle gate in the four turbulent conditions and without turbulence.</sub></p>

### 3. Online estimate on real imagery

<p align="center"><img src="docs/figures/real.png" width="100%"></p>
<p align="center"><sub><b>Figure 4.</b> Tilt variance estimated by the unchanged estimator against the scintillometer reading, 119 one-minute records over a 724 m path (data of Saha et al., not redistributed here).</sub></p>

| Criterion (declared before processing) | Value [95% CI] | Outcome |
|:---|---:|:---|
| Rank correlation with scintillometer | 0.863 [0.585, 0.903] | passed |
| Log-log slope | 0.833 [0.665, 0.883] | passed |
| Median ratio, estimated / measured Cn2 | 0.399 [0.366, 0.413] | passed |
| Median R2 of the structure-function fit | 0.615 | failed |

### 4. Ablation: tilt against blur

Share of the ATE increase of the baseline that is reproduced when only one component of the degradation is
applied (45 units). The shares sum to more than one at D/r0 = 16; the effects are sub-additive.

| Level | Share of Tilt Alone | Share of Blur Alone |
|:---|---:|---:|
| D/r0 = 8 | 0.818 | 0.087 |
| D/r0 = 16 | 0.836 | 0.362 |

## Protocol

| Item | Setting |
|:---|:---|
| Data | EuRoC MH_01 to MH_05, V1_01, V1_02, V2_01, V2_02; up to three windows of 400 frames per sequence |
| Unit | one (sequence, seed) pair; ATE averaged over the windows of a unit |
| Metric | ATE of the causal trajectory after rigid alignment over the window; position ANEES |
| Turbulence | D = 0.10 m, 550 nm, 5.86 µrad per pixel; quasi-static (coherence time = frame interval) and fast (5 m/s) |
| Statistics | paired median ratio, 10 000 bootstrap resamples, Wilcoxon signed-rank, Holm adjustment, TOST at ±15% |
| Pre-registration | `preregistration/ON_KAYIT_FAZ2.md` with seven dated deviations and a SHA-256 history |

## Installation

```bash
# creates the conda environment `pyeco` (Python 3.12, torch 2.5.1 + CUDA 12.4, gtsam 4.2.1, numpy 1.26.4)
bash install.sh
# machine without a GPU
bash install.sh --cpu
```

Reference machine: Ubuntu 24.04, RTX 3090. The numpy pin matters: the gtsam wheel fails with numpy 2.

## Data

```bash
# EuRoC sequences in ASL format under one root: $SLAM_DATA_ROOT/MH01/mav0/...
export SLAM_DATA_ROOT=/data/sdata
```

EuRoC MAV is distributed by ETH Zurich (Burri et al., 2016) and is not included; Figure 1 shows one frame of
MH_05 for illustration. The real imagery of Saha et al. (2022) and the AirSim scene are not redistributed;
only derived per-minute records and analysis outputs are stored. The learned restoration networks used in
one experiment are not included (no code, no weights).

## Running

```bash
# plan of an experiment without running it
python sapma7_kos.py --kuru
# adaptive-gate confirmatory grid, one part (resumable)
python kapi_dogrula.py --parca II --dokumler DUMP_DIR
# main grid on one GPU: dry check, then start with 10 workers
./faz2_kos.sh kuru
./faz2_kos.sh baslat 10
```

The backend options of the arms are ordinary command-line flags, for example the adaptive gate:
`--turb_cal_factor 1.0 --turb_dimm_sigma_scale 1.0 --turb_gate 3` on the arm `A0C+A3`.
Identifiers, comments, and run logs of the experiment scripts are in Turkish (`kos` = run, `analiz` =
analysis, `kapi` = gate, `hizli` = fast regime, `onarim` = restoration, `sapma` = deviation, `ontest` =
pre-test, `kuru` = dry run).

## Reproducing the paper

```bash
# code fingerprints of both trees
python scripts/verify_identity.py
# headline numbers of the paper against the analysis files
python scripts/check_paper_numbers.py
# tables of this README against the record files
python scripts/readme_tables.py --check
# figures of this README
python scripts/make_readme_figures.py
# decision tables from the stored records (each prints its hypotheses and outcomes)
python -m turb.faz2_analiz
python -m turb.h2_analiz
python hizli_analiz.py
python a3s_analiz.py
python onarim_analiz.py
python kapi_dogrula_analiz.py
```

Every record in `results/` carries the fingerprint of the code that produced it (`docs/IDENTITY.md`).

## Repository layout

| Path | Content |
|:---|:---|
| `turb/` | simulator, online estimator, covariance arms, metrics, runners, analyses |
| `gtsam_backend.py`, `stereo_front_end.py`, `realtime_pipeline.py` | the odometry pipeline (PyEco-SLAM) |
| `kapi_*.py`, `sapma7_*.py`, `hizli_*.py`, `a3s_*.py`, `onarim_*.py` | runners and analyses of the dated deviations |
| `airsim/` | recording and grid scripts of the rendered test bed |
| `results/` | records and analysis outputs of every experiment (`evidence` is a link to it) |
| `preregistration/` | pre-registration document, earlier versions, deviations, hash history, pre-test reports |
| `reference/main_tree/` | the three files that differ in the main tree, and the patch |
| `scripts/` | identity check, number checks, README tables and figures |
| `docs/` | identity, provenance, figures |

## License

MIT for the code of this repository (`LICENSE`). Vendored third-party code keeps its own license (`NOTICE`).

## Citation

```bibtex
@unpublished{asil2026turbvo,
  author = {Asil, Ufuk},
  title  = {Stereo Visual Odometry Under Atmospheric Optical Turbulence:
            Online Image-Motion Estimation and Adaptive Outlier Gating},
  note   = {Manuscript},
  year   = {2026}
}
```
