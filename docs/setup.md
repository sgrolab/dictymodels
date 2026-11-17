## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/sgrolab/dictymodels.git
cd dictymodels
```

### 2. Install uv
The project's dependencies are pinned in `pyproject.toml` and `uv.lock` and managed with
[uv](https://docs.astral.sh/uv/). If you don't already have it:

**macOS/Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows:**
```powershell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternative (using pip):**
```bash
pip install uv
```

### 3. Create the Environment
```bash
uv sync
```
This creates a virtual environment in `.venv` with Python 3.12 (downloading it if necessary)
and installs the pinned dependencies from `uv.lock`. They are the versions that were current
in August 2024, when the model and figure code was last changed.

There is no need to activate the environment: prefix every command with `uv run`, which runs
it inside `.venv`:
```bash
uv run python --version
```

If you prefer another environment manager (for example conda), install the dependencies
listed in `pyproject.toml` and drop the `uv run` prefix from the commands below.

### 4. Point the Code at Your Data Directory
The simulations write their results outside the repository, in a directory the code calls
`PROJECT_DIR`. It is read from the `ROOT_DIR` environment variable by
`compare_models/config.py`, which also loads a `.env` file from the repository root if one
exists (`.env` is git-ignored). Either create that file once:
```bash
echo "ROOT_DIR=/path/to/your/data/directory" > .env
```
or set the variable for each command:
```bash
ROOT_DIR=/path/to/your/data/directory uv run compare_models/plotFigures.py
```
On Windows, copy `.env.example` to `.env` and edit it in a text editor instead.

Give the full path: a relative one is taken from whichever directory a command is run in.
The directory and the folders inside it are created when they are first needed. Each script
writes its results there as Python pickle files, in a folder of its own.

### 5. Verify the Installation
```bash
uv run compare_models/SingleCell_singleSpike.py
uv run compare_models/plotFigures.py --save figs figure_02
```
The first command takes about a second and writes five files to
`<ROOT_DIR>/sc_singleSpike_data/`. The second draws Figure 2 of the paper from them and
saves it as `figs/figure_02.png`. If either stops with `RuntimeError: ROOT_DIR is not set`,
go back to step 4.

## Running the Code

### Simulations
Each experiment is one script in `compare_models/`. It runs the models, with the parameters
written out at the top of each model's section, and saves the results for `plotFigures.py`
to read. A script is run directly:
```bash
uv run compare_models/SingleCell_adaptiveSpike.py
```
`compare_models/runSimulations.py` runs them in bulk, each in its own process:
```bash
uv run compare_models/runSimulations.py --dry-run     # list what a full run would do
uv run compare_models/runSimulations.py --jobs 4      # run everything, four at a time
uv run compare_models/runSimulations.py SingleCell    # the scripts whose names start so
uv run compare_models/runSimulations.py SingleCell_singleSpike.py Population_add_cAMP_SCnoise.py
uv run compare_models/runSimulations.py --help
```
A full run includes the sweep described below, which is nearly all of its work; `--jobs`
shares the work among processes that run at the same time. Scripts that fail are reported
at the end and the rest still run, unless `--stop-on-error` is given. `--log` keeps each
script's output in a file.

| Figure | Script                              | Writes to                      | Time    |
|--------|-------------------------------------|--------------------------------|---------|
| 2      | `SingleCell_singleSpike.py`         | `sc_singleSpike_data/`         | 1 s     |
| 3      | `Population_add_cAMP_SCnoise.py`    | `popAddcAMPscNoise_data/`      | 1 min   |
| 4, S1  | the sweep, described below          | `pop_firingRate_data/`         | 40 h    |
| 5      | `SingleCell_adaptiveSpike.py`       | `sc_adaptiveSpike_data/`       | 1 s     |
| 5      | `SingleCell_bifurcationDynamics.py` | `sc_bifurcationDynamics_data/` | 1 s     |
| 6, S2  | `SingleCell_FoldChangeDetection.py` | `sc_foldChange_data/`          | 6.5 min |
| S3     | `SingleCell_stepVsRamp.py`          | `sc_stepVsRamp_data/`          | 3 s     |
| S4     | `SingleCell_Entrainment.py`         | `sc_Entrainment_data/`         | 1 min   |
| S5     | `SingleCell_timescaleSeparation.py` | `sc_timescaleSeparation_data/` | 1 s     |

The times are from one run on a four-core Linux machine. `Population_add_cAMP_SCnoise.py`
needs about 7 GB of memory and writes 5.3 GB; the other scripts in the table need less than
1 GB and write 65 MB between them. `Population_add_cAMP.py` is the Figure 3 experiment
without noise in the single cells; no figure uses its results, so a full run leaves it out.

Noise is drawn from a random number generator with a fixed seed, so running a script again
reproduces its results.

The scripts are divided into cells by `#%%` markers, so they can also be run one cell at a
time in Spyder or VS Code. Set the interpreter to `.venv/bin/python`
(`.venv\Scripts\python.exe` on Windows) and the working directory to `compare_models/`, so
that the scripts can import the models and `config.py`.

### The Density and Flow-Rate Sweep
Figures 4 and S1 show how often a population of 100 cells fires at each point of a grid of
cell densities and flow rates, for every model, with and without noise in the single cells.
There are twelve scripts: `Population_DensityFlowSweep_<Model>_Noise.py` and
`..._noNoise.py` for each of the five models, and two more for Kamino ending in `_log`, in
which the flow rates are spaced logarithmically. Each run simulates one grid point, chosen
by two indices on the command line, and writes one file:
```bash
uv run compare_models/Population_DensityFlowSweep_Sgro_Noise.py 13 20
# writes <ROOT_DIR>/pop_firingRate_data/noise/Sgro/j_13_rho_20.pickle
```

| Model     | Scripts | Grid    | Indices | Each run      | All runs       |
|-----------|---------|---------|---------|---------------|----------------|
| Sgro      | 2       | 26 x 26 | 0 to 25 | 2 s, 0.1 GB   | 41 min, 2.3 GB |
| Gregor    | 2       | 26 x 26 | 0 to 25 | 2 s, 0.1 GB   | 37 min, 2.6 GB |
| Goldbeter | 2       | 25 x 25 | 0 to 24 | 4 s, 0.1 GB   | 1.4 h, 3.3 GB  |
| Maeda     | 2       | 26 x 26 | 0 to 25 | 37 s, 0.2 GB  | 14 h, 23 GB    |
| Kamino    | 4       | 31 x 31 | 0 to 30 | 21 s, 0.2 GB  | 23 h, 77 GB    |

"Each run" is the time and peak memory of one grid point. "All runs" is the time on one
processor core and the disk space for all of that model's scripts. In all that is 9,150
runs, about 40 hours and 110 GB. One grid point of each script was timed, so the totals are
estimates.

`runSimulations.py` runs a sweep script once for every point of its grid. When all twelve
scripts have a file for every grid point, it runs
`Population_DensityFlowSweep_Consolidate.py`, which collects the firing rates into one small
file per script for the figures:
```bash
uv run compare_models/runSimulations.py Population_DensityFlowSweep --jobs 4
```
Ctrl-C stops it, and running it again with `--resume` carries on with the grid points that
have no file yet. A grid point that fails is named in the summary and its file removed, so
`--resume` tries it again. On a four-core machine, three processes ran the 676 grid points
of one script in under 8 minutes, against about 21 minutes for one process.

The lab ran the sweep on an LSF cluster instead. `Population_DensityFlowSweep_GenerateBsubs.py`
writes the job lists for that: twelve text files in `compare_models/`, one per script, with
one `bsub` command per grid point. It adds to files that already exist, so delete old lists
before running it again. After a sweep run that way, run the consolidate script by hand:
```bash
uv run compare_models/Population_DensityFlowSweep_Consolidate.py
```
It reads all twelve result folders and stops if one is missing. A grid point that has no
file is left at zero without a warning, so first check that each folder holds one file per
grid point (676, 625 or 961). `runSimulations.py` makes that check itself.

### Figures
`compare_models/plotFigures.py` draws Figures 2 to 6 and S1 to S5 of the paper from the data
in `ROOT_DIR`. Without arguments it opens each figure in a window; with `--save` it writes
them to disk instead:
```bash
uv run compare_models/plotFigures.py                          # show every figure
uv run compare_models/plotFigures.py figure_03 figure_S04     # just these
uv run compare_models/plotFigures.py --save figs --format pdf --dpi 300
uv run compare_models/plotFigures.py --help
```
A figure whose data has not been generated yet is reported as failed and the others are
still drawn. Figure 3 loads 5.3 GB of traces and needs about 7 GB of memory.

Each figure is a function, `figure_02` ... `figure_S05`. In Spyder or VS Code, with the
interpreter and working directory set as described under Simulations, `import plotFigures`
and then call one, e.g. `plotFigures.figure_03()`.

### Using a Model in Your Own Script
Each `*Func.py` file defines a class for a single cell and classes for a population. A
script saved in `compare_models/` can import them:
```python
import numpy as np
import matplotlib.pyplot as plt
import Kamino2017Func as kf

dt = 0.0001                       # time step
t = np.arange(0, 30, dt)          # time axis, in the model's own units
signal = np.zeros(len(t))         # external cAMP: a step from 0 to 1 at t = 5
signal[len(t) // 6:] = 1

params = {'tau': 1.5, 'n': 2, 'K': 4, 'kt': 2, 'delta': 0.01, 'gamma': 3, 'rho': 0.01}
cell = kf.Cell([0.01, 0.06, 0.005], params, t)    # initial x, y and z
cell.run(dt, signal)

plt.plot(cell.t, cell.y)          # y is the variable the figures show for this model
plt.xlabel('time')
plt.ylabel('y')
plt.show()
```
The arguments differ from model to model; `SingleCell_singleSpike.py` runs each of the five
once and is the shortest example of them all. To put the models on common scales, the
experiment scripts divide time by `Nt_<model>` and the response by `Nh_<model>`, the
normalization constants that every script reads at its top. Figure 2 of the paper shows how
they were chosen.

The classes for a population of noisy cells (`Population_scNoise`, `Gregor2010_pop` and
`Sgro2015_pop`) keep only the mean of each variable over the cells, for example
`pop.y_mean`, because every variable of every cell at every time step can take gigabytes.
To look inside a run, call `run(..., keep_cells=True)`: each variable is then also kept for
every cell, as `pop.y`, an array of cells by time steps. `Population_add_cAMP_SCnoise.py`
does this for Figure 3, which shows single cells. The noise can be given to `run()` whole,
as an array of time steps by cells, or as a function that returns the next rows of one, so
that it is drawn while the run proceeds, as the sweep scripts do. See `popHistory.py`.
