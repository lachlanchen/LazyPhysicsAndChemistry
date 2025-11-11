# Agent Guide for QUANTUM

This repository contains helpers and inputs for running Gaussian jobs and related tooling. Please follow these conventions when making changes or automating tasks here.

## Repo Structure Notes

- `Gaussian/` is a symlink pointing outside this folder (to `../Gaussian/`). Any changes under `Gaussian/…` modify files in the target directory, not just this repo. Be deliberate when editing under this path.
- `Gaussian/run_gaussian.sh` runs Gaussian 16 (`g16`) on a given `.com`/`.gjf` input and optionally opens the result in GaussView.
- `Gaussian/_test_run/` contains a tiny water example used to sanity‑check the install.
- `multiwfn/` holds Multiwfn source/docs; unrelated to Gaussian runner changes unless explicitly working on Multiwfn.

## Coding and Editing Conventions

- Keep changes minimal and focused on the task at hand.
- Do not add heavy dependencies or networked installers.
- Prefer small, self‑contained shell scripts for automation.
- When touching anything in `Gaussian/`, remember it crosses the symlink boundary; avoid accidental reorganizations there.

## Gaussian Runner Usage

- Run a job and open GaussView:
  - `Gaussian/run_gaussian.sh <path/to/input.com>`
- Suppress GaussView:
  - `Gaussian/run_gaussian.sh --no-view <path/to/input.com>`
- Specify binaries explicitly:
  - `Gaussian/run_gaussian.sh --g16 ~/gaussian/g16/g16 --gview ~/gaussian/gv/gview_safe.sh <input>`

Behavior:
- Writes `<basename>.log` next to the input.
- Detects `%chk=…` in the input; if the checkpoint exists, GaussView opens the `.chk`, otherwise the `.log`.
- Uses `GAUSS_SCRDIR` if set, else defaults to `~/gaussian/scr`.

## Environment Assumptions

- Gaussian 16 installed at `~/gaussian/g16/g16`, or `g16` on `PATH`.
- Scratch directory set via `GAUSS_SCRDIR` (the runner will create `~/gaussian/scr` if missing).
- GaussView installed under `~/gaussian/gv/`. The runner prefers `gview_safe.sh` if present, then `gview.sh`.

Recommended GaussView wrapper (stabilizes GL/Wayland):

```
#!/usr/bin/env bash
set -euo pipefail
GV_SH="$HOME/gaussian/gv/gview.sh"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-xcb}"
export LIBGL_ALWAYS_SOFTWARE="${LIBGL_ALWAYS_SOFTWARE:-1}"
export __GLX_VENDOR_LIBRARY_NAME="${__GLX_VENDOR_LIBRARY_NAME:-mesa}"
exec "$GV_SH" "$@"
```

Save as `~/gaussian/gv/gview_safe.sh` and make it executable.

## Troubleshooting

- Success criterion: `Normal termination of Gaussian` appears near the end of the `.log`.
- If GaussView fails to launch in Wayland/remote sessions, use the safe wrapper above and pass `--gview ~/gaussian/gv/gview_safe.sh` to the runner.
- If scratch space errors occur, verify free disk under `GAUSS_SCRDIR` and permissions.

