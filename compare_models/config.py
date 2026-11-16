"""
Locations of the code, the experimental data and the simulation data.

Simulation data live outside the repository, in the directory named by the ROOT_DIR
environment variable. Set it in a file named .env in the repository root (copy
.env.example; .env is git-ignored), or in the shell.
"""
import os
from pathlib import Path
from dotenv import load_dotenv

# model functions and scripts, and the experimental data, are found from this file
CODE_DIR = Path(__file__).resolve().parent
EXP_DATA_DIR = CODE_DIR.parent / 'exp_data'

# simulation data are found from ROOT_DIR
load_dotenv(CODE_DIR.parent / '.env')

if not os.getenv('ROOT_DIR'):
    raise RuntimeError('ROOT_DIR is not set: copy .env.example to .env in the repository '
                       'root and put the path of your data directory in it')

PROJECT_DIR = Path(os.getenv('ROOT_DIR')).expanduser().resolve()

def output_dir(name):
    """Return the directory PROJECT_DIR/name, creating it if it does not exist."""
    path = PROJECT_DIR / name
    path.mkdir(parents=True, exist_ok=True)
    return path
