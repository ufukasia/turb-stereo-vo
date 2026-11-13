# Code identity

Every result record in `results/` carries `code_fp`, the first 12 hex digits of a SHA-1 over the
names and SHA-1 digests of the source files listed in `turb/runner.py:code_fingerprint()`.

| Tree | `code_fp` | Used by |
|---|---|---|
| main tree | `08e28b83407c` | main grid, factorial, fast regime, scalar variant, restoration, probes |
| gate tree | `bedca2849bb2` | adaptive gate (deviation 6) and gated generic baseline (deviation 7) |

The repository root is the gate tree. It is a superset of the main tree: with `--turb_gate 0`
(default) the gate code is inactive. Three fingerprinted files differ between the trees:
`config.py`, `gtsam_backend.py`, `turb/backend.py`. The main-tree versions are kept byte for byte in
`reference/main_tree/`, and the difference is in `reference/main_tree_to_gate_tree.patch`.

```bash
# recompute both fingerprints (no dependencies beyond the standard library)
python scripts/verify_identity.py
```

Expected output: `gate tree: bedca2849bb2 ... OK` and `main tree: 08e28b83407c ... OK`.
To run main-tree experiments exactly, copy `reference/main_tree/*` over the root files.

## Fingerprinted files are frozen

The fingerprinted files must not be edited, reformatted, or re-encoded; any change alters `code_fp`
and breaks the link to the stored records. For this reason one machine-specific default path was left
in place in a fingerprinted file:

- `turb/runner.py`, line 219: default of `TURB_DUMP_CACHE` (override with the environment variable).

In all other files local paths were replaced by generic defaults (`/data/...`); set `SLAM_DATA_ROOT`
and the command-line options instead. Result records and pre-test reports are stored as written and
may contain paths of the machine on which they were produced.

## Pre-registration

`preregistration/ON_KAYIT_FAZ2.md` (Turkish) has SHA-256
`47c481c7a23983f473268b2f6a65c5e2ae9feeedf02a5a33e4007ef0c4a4aeb8`; the hash history after each dated
deviation is in `preregistration/on_kayit_sha256.txt`.

```bash
sha256sum preregistration/ON_KAYIT_FAZ2.md
```
