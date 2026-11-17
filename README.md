Robust coordination of collective oscillatory signaling requires single-cell excitability
and fold-change detection
=========================================================================================

## Project Abstract
Complex multicellular behaviors are coordinated at the level of biochemical signaling
networks, yet how this decentralized mechanism enables robust control in variable
environments and over many orders of magnitude of spatiotemporal scales remains an open
question. A stunning example of these behaviors is found in the microbe *Dictyostelium
discoideum*, which uses the small molecule cyclic AMP (cAMP) to drive the propagation of
collective signaling oscillations leading to multicellular development. The critical design
features of the *Dictyostelium* signaling network remain unclear despite decades of
mathematical modeling and experimental interrogation because each model makes different
assumptions about the network architecture and in general, normalizing models for direct
comparison presents a major challenge. We overcome this challenge by using recent
experimental data to normalize the time and response scales of five major signal relay
network models to one another and assess their ability to recapitulate
experimentally-observed population and single-cell dynamics. We find that to successfully
reproduce the full range of observed dynamical behaviors, single cells must be excitable and
respond to the relative fold-change of environmental signals. This suggests these features
represent robust principles for coordinating cellular populations through oscillatory
signaling and that single-cell excitable dynamics are a generalizable route for controlling
population behaviors.

## Models
This repository holds the five models, the simulations that compare them and the script
that draws the figures of the paper.

| Label in the figures     | Model of                 | File                     |
|--------------------------|--------------------------|--------------------------|
| Receptor Desensitization | Martiel & Goldbeter 1987 | `Goldbeter1987Func.py`   |
| CDINFB                   | Maeda et al. 2004        | `MaedaLoomis2004Func.py` |
| Phase Oscillator         | Gregor et al. 2010       | `Gregor2010Func.py`      |
| IPNFB                    | Sgro et al. 2015         | `Sgro2015Func.py`        |
| IFFL                     | Kamino et al. 2017       | `Kamino2017Func.py`      |

CDINFB stands for coupled direct and indirect negative feedback, IPNFB for interlocking
positive-negative feedback, and IFFL for incoherent feedforward loop.

## Repository Structure
- `docs/` - Installation and usage instructions (`setup.md`)
- `compare_models/*Func.py` - The five models, one file each, with classes for a single cell
  and for a population of cells
- `compare_models/popHistory.py` - What a run of a population keeps (the mean over the cells
  by default, every cell on request), how it is given its noise, what the loops that numba
  compiles share, and a fast way to find the peaks of a sweep run
- `compare_models/SingleCell_*.py`, `compare_models/Population_add_cAMP*.py` - One script per
  experiment; each writes its results to the data directory
- `compare_models/Population_DensityFlowSweep_*.py` - The sweep over cell density and flow
  rate behind Figures 4 and S1, one grid point per run
- `compare_models/runSimulations.py` - Runs those scripts in bulk, as many at a time as the
  computer has processors, and the sweep scripts once for every grid point (`--help` for
  options)
- `compare_models/plotFigures.py` - Draws Figures 2 to 6 and S1 to S5 from those results
  (`--help` for options)
- `compare_models/config.py` - Tells the scripts where the experimental data and the data
  directory are
- `exp_data/` - Published experimental measurements that the figures show beside the
  simulations
- `Notebooks/` - Tutorials written for an earlier version of the code; they do not run with
  the current one

## Requirements
- Python 3.12
- The Python packages pinned in `pyproject.toml` and `uv.lock` (numpy, scipy, matplotlib,
  numba, pandas, openpyxl, python-dotenv), installed with `uv sync` and used through
  `uv run`
- A data directory, pointed to by the `ROOT_DIR` environment variable, where the simulations
  write their results
- About 3 GB of free memory to draw Figure 3, and 0.6 GB of disk for its data; the other
  figures outside the sweep need less
- For Figures 4 and S1, the sweep: 9,150 runs that add up to about 4 hours on one processor
  core (`runSimulations.py` runs one on every core at a time), with about 0.5 GB of memory
  for each core and about 110 GB of disk (see `docs/setup.md`)

## Installation and Usage
See `docs/setup.md` for installing the environment, configuring the data directory, running
the simulations and plotting the figures. In short:
```bash
uv sync
echo "ROOT_DIR=/path/to/your/data/directory" > .env
uv run compare_models/runSimulations.py --dry-run     # list the simulations
uv run compare_models/runSimulations.py SingleCell    # run the single-cell ones
uv run compare_models/plotFigures.py --save figs      # draw the figures that have data
```

## Citation
Huyan C, Golden A, Zhu X, Mehta P, Sgro AE. Robust coordination of collective oscillatory
signaling requires single-cell excitability and fold-change detection. bioRxiv (2021).
https://doi.org/10.1101/2021.09.02.457527

## Contact
Contact Allyson Sgro at sgroa@janelia.hhmi.org.
