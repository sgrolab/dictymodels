#!/usr/bin/env python
"""Run the simulation scripts that generate the data for the figures.

Each script is run in its own Python process, exactly as it would be by hand. As many run
at a time as the computer has processors, or as --jobs says. The scripts of the density and
flow-rate sweep simulate one grid point per run, so each of them is run once for every point
of its grid; that is 9,150 runs and nearly all of the work. When every grid point of the
sweep has its file, the consolidate script collects them for the figures.

Without arguments every script is run except those in EXCLUDE; name scripts, or the start
of the names of several, to run a subset. The data directory is the one the scripts use
(ROOT_DIR from the environment or .env).

Ctrl-C stops a run; --resume then carries on with the grid points that have no file yet.
"""
import argparse
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
SWEEP_PREFIX = 'Population_DensityFlowSweep_'
CONSOLIDATE = SWEEP_PREFIX + 'Consolidate.py'

# Scripts left out of a full run, and why; they can still be named on the command line
EXCLUDE = {'Population_add_cAMP.py': 'no figure uses its results'}

# Scripts with the names of simulations that are not simulations
NOT_SIMULATIONS = {SWEEP_PREFIX + 'GenerateBsubs.py'}

# The sweep scripts, by their names after SWEEP_PREFIX. Each takes two indices on the command
# line, simulates that point of its grid and writes one file under pop_firingRate_data/.
# Here are the number of values of each index and the file of a grid point, which repeat
# what the scripts themselves say: a run whose file does not appear is reported as failed.
SWEEP = {
    'Sgro_Noise':         (26, 'noise/Sgro/j_%.2i_rho_%.2i.pickle'),
    'Sgro_noNoise':       (26, 'no_noise/Sgro/j_%.2i_rho_%.2i.pickle'),
    'Gregor_Noise':       (26, 'noise/Gregor/rho_%.2i_k_%.2i.pickle'),
    'Gregor_noNoise':     (26, 'no_noise/Gregor/rho_%.2i_k_%.2i.pickle'),
    'Goldbeter_Noise':    (25, 'noise/Goldbeter/kc_%.2i_h_%.2i.pickle'),
    'Goldbeter_noNoise':  (25, 'no_noise/Goldbeter/kc_%.2i_h_%.2i.pickle'),
    'Maeda_Noise':        (26, 'noise/Maeda/g_%.2i_rho_%.2i.pickle'),
    'Maeda_noNoise':      (26, 'no_noise/Maeda/g_%.2i_rho_%.2i.pickle'),
    'Kamino_Noise':       (31, 'noise/Kamino/j_%.2i_rho_%.2i.pickle'),
    'Kamino_noNoise':     (31, 'no_noise/Kamino/j_%.2i_rho_%.2i.pickle'),
    'Kamino_Noise_log':   (31, 'noise/Kamino_log/j_%.2i_rho_%.2i.pickle'),
    'Kamino_noNoise_log': (31, 'no_noise/Kamino_log/j_%.2i_rho_%.2i.pickle'),
}


def all_scripts():
    """The simulation scripts in the order they are run: single cells, populations, the
    sweep, and last the script that consolidates the sweep."""
    def named(pattern):
        return sorted((p.name for p in CODE_DIR.glob(pattern)), key=str.lower)
    sweep = [SWEEP_PREFIX + name + '.py' for name in SWEEP]
    found = set(named(SWEEP_PREFIX + '*.py')) - NOT_SIMULATIONS - {CONSOLIDATE}
    if found != set(sweep):
        sys.exit('SWEEP does not match the sweep scripts: '
                 + ' '.join(sorted(found ^ set(sweep))))
    populations = [s for s in named('Population_*.py') if not s.startswith(SWEEP_PREFIX)]
    return named('SingleCell_*.py') + populations + sweep + [CONSOLIDATE]


def resolve(targets):
    """The scripts that the targets name, in running order. A target is a script, with or
    without .py or a folder before it, or the start of the names of several scripts."""
    scripts = all_scripts()
    if not targets:
        return [s for s in scripts if s not in EXCLUDE]
    chosen = set()
    for target in targets:
        stem = Path(target).name.removesuffix('.py')
        matches = [s for s in scripts if s == stem + '.py'] or \
                  [s for s in scripts if s.startswith(stem)]
        if not matches:
            sys.exit(f'unknown target: {target}\nthe scripts are:\n  ' + '\n  '.join(scripts))
        chosen.update(matches)
    return [s for s in scripts if s in chosen]


def grid_points(script, data_dir):
    """The runs of a sweep script: for each, its two indices and the file it writes."""
    points, name = SWEEP[script[len(SWEEP_PREFIX):-3]]
    return [((i, j), data_dir / 'pop_firingRate_data' / (name % (i, j)))
            for i in range(points) for j in range(points)]


def processors():
    """How many processors this process may use."""
    try:
        return len(os.sched_getaffinity(0))     # Linux: what taskset or a cluster job allows
    except AttributeError:                      # macOS and Windows cannot say
        return os.cpu_count() or 1


def hms(seconds):
    seconds = int(round(seconds))
    return f'{seconds // 3600}:{seconds // 60 % 60:02d}:{seconds % 60:02d}'


class Runner:
    """Runs scripts in processes of their own, and stops them when told to."""

    def __init__(self, log_dir, stop_on_error):
        self.log_dir = log_dir
        self.stop_on_error = stop_on_error
        self.stop = threading.Event()   # set to start no more runs
        self.lock = threading.Lock()
        self.processes = set()

    def run(self, script, indices=(), output=None):
        """Run a script, with the indices of a grid point if it is a sweep script. Returns
        a description of how it went; 'code' is None if it was not started."""
        result = {'script': script, 'indices': indices, 'code': None, 'seconds': 0.0,
                  'text': ''}
        start = time.monotonic()
        with self.lock:
            if self.stop.is_set():
                return result
            process = subprocess.Popen(
                [sys.executable, str(CODE_DIR / script)] + [str(i) for i in indices],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace')
            self.processes.add(process)
        result['text'] = process.communicate()[0]
        with self.lock:
            self.processes.discard(process)
        result['code'] = process.returncode
        result['seconds'] = time.monotonic() - start
        if output is not None:
            if result['code'] != 0:
                try:
                    output.unlink(missing_ok=True)  # whatever it wrote is not to be trusted
                except OSError as error:
                    result['text'] += f'could not remove {output}: {error}\n'
            elif not output.is_file():
                result['code'] = 'no file'
                result['text'] += f'the run ended normally but did not write {output}\n'
        if self.stop_on_error and result['code'] != 0:
            self.stop.set()     # here, before another run can be started
        if self.log_dir is not None and (result['text'] or not indices):
            name = '_'.join([script[:-3]] + [str(i) for i in indices]) + '.log'
            (self.log_dir / name).write_text(result['text'])
        return result

    def interrupt(self):
        """Start no more runs and end those that are under way."""
        self.stop.set()
        with self.lock:
            for process in self.processes:
                process.terminate()


def show_output(text, lines=20):
    """Print the end of what a run wrote, indented."""
    for line in text.rstrip().split('\n')[-lines:] if text.strip() else []:
        print('    ' + line)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('targets', nargs='*', metavar='TARGET',
                        help='a script (SingleCell_singleSpike.py) or the start of the names '
                             'of several (SingleCell, Population_DensityFlowSweep_Kamino); '
                             'default: every script except ' + ' '.join(EXCLUDE))
    parser.add_argument('--jobs', metavar='N', type=int, default=processors(),
                        help='how many processes to run at a time (default: one for each '
                             'processor, %(default)s here)')
    parser.add_argument('--resume', action='store_true',
                        help='skip the grid points of the sweep whose file is already there')
    parser.add_argument('--dry-run', action='store_true',
                        help='list what would run, in order, and exit')
    parser.add_argument('--stop-on-error', action='store_true',
                        help='start nothing more after a failure (default: carry on and '
                             'report)')
    parser.add_argument('--log', metavar='DIR', type=Path,
                        help="write each script's output to DIR/<script>.log, and that of a "
                             'sweep run that printed anything to DIR/<script>_<i>_<j>.log')
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error('--jobs must be at least 1')

    scripts = resolve(args.targets)
    try:
        from config import PROJECT_DIR
    except RuntimeError as error:   # ROOT_DIR is not set
        parser.error(str(error))

    # what there is to do: scripts that run once, and the grid points of the sweep scripts
    once = [s for s in scripts if not s.startswith(SWEEP_PREFIX)]
    sweep = {s: grid_points(s, PROJECT_DIR) for s in scripts
             if s.startswith(SWEEP_PREFIX) and s != CONSOLIDATE}
    todo = {s: [(indices, output) for indices, output in points
                if not (args.resume and output.is_file())] for s, points in sweep.items()}
    total = sum(len(points) for points in todo.values())
    jobs = max(1, min(args.jobs, len(once) + total))    # no more processes than runs

    print(f'data directory: {PROJECT_DIR}')
    print(f'{len(scripts)} script(s), {jobs} process(es) at a time:')
    for s in scripts:
        note = ''
        if s in sweep:
            points = SWEEP[s[len(SWEEP_PREFIX):-3]][0]
            there = len(sweep[s]) - len(todo[s])
            note = f'{len(todo[s])} runs ({points} x {points}' + \
                   (f', {there} already there)' if args.resume else ')')
        elif s == CONSOLIDATE:
            note = 'once every grid point of the sweep has its file'
        print(f'  {s:<51}{note}'.rstrip())
    if not args.targets:
        for s, why in EXCLUDE.items():
            print(f'left out: {s} ({why})')
    if sweep:
        print(f'{total:,} sweep run(s) in all')
    if args.dry_run:
        return 0
    if args.log is not None:
        args.log.mkdir(parents=True, exist_ok=True)

    runner = Runner(args.log, args.stop_on_error)
    began = time.monotonic()
    reported = began
    interrupted = False
    print()
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(runner.run, s) for s in once]
        futures += [pool.submit(runner.run, s, indices, output)
                    for s, points in todo.items() for indices, output in points]
        try:
            done = failed = 0
            for future in as_completed(futures):
                result = future.result()
                name = ' '.join([result['script']] + [str(i) for i in result['indices']])
                if result['code'] is None:
                    continue
                if result['code'] != 0:
                    print(f'FAILED  {hms(result["seconds"])}  {name}  ({result["code"]})')
                    show_output(result['text'])
                elif not result['indices']:
                    print(f'ok      {hms(result["seconds"])}  {name}')
                    show_output(result['text'])
                if result['indices']:
                    done += 1
                    failed += result['code'] != 0
                    if done == total or time.monotonic() - reported > 60:
                        reported = time.monotonic()
                        print(f'  {done:,} of {total:,} sweep runs done, {failed:,} failed, '
                              f'{hms(reported - began)} elapsed', flush=True)
        except KeyboardInterrupt:
            interrupted = True
            print('\ninterrupted: stopping the runs under way')
            runner.interrupt()
    results = [future.result() for future in futures]   # all have ended by now

    # the consolidate script, if it was asked for and the sweep is all there
    consolidated = None
    if CONSOLIDATE in scripts and not runner.stop.is_set():
        missing = {s: sum(not output.is_file() for _, output in grid_points(s, PROJECT_DIR))
                   for s in all_scripts() if s.startswith(SWEEP_PREFIX) and s != CONSOLIDATE}
        if any(missing.values()):
            consolidated = 'not run: no file for ' + ', '.join(
                f'{n:,} grid point(s) of {s}' for s, n in missing.items() if n)
        else:
            consolidated = runner.run(CONSOLIDATE)

    print("\n===== summary (the time of a sweep script is that of its runs added up) =====")
    bad = 0
    for s in scripts:
        mine = [r for r in results if r['script'] == s]
        ran = [r for r in mine if r['code'] is not None]
        failures = [r for r in ran if r['code'] != 0]
        seconds = hms(sum(r['seconds'] for r in ran))
        if s == CONSOLIDATE:
            if isinstance(consolidated, dict):
                ok = consolidated['code'] == 0
                print(f'{"ok    " if ok else "FAILED"}  {hms(consolidated["seconds"])}  {s}')
                if not ok:
                    show_output(consolidated['text'])
            else:
                ok = False
                print(f'not run          {s}')
                if consolidated:
                    print(f'    {consolidated}')
        elif s in sweep:
            ok = not failures and len(ran) == len(mine)
            line = f'{len(ran) - len(failures):,} of {len(todo[s]):,} runs ok'
            if failures:
                line += '; failed: ' + ', '.join(
                    ' '.join(str(i) for i in r['indices']) for r in failures[:10])
                line += ', ...' if len(failures) > 10 else ''
            state = 'ok    ' if ok else 'FAILED' if failures else 'not run'
            print(f'{state}  {seconds}  {s}  {line}')
        else:
            ok = bool(ran) and not failures
            print(f'{"ok    " if ok else "FAILED" if ran else "not run"}  {seconds}  {s}')
        bad += not ok
    print(f'{len(scripts) - bad} script(s) ok, {bad} not; '
          f'{hms(time.monotonic() - began)} elapsed')
    return 130 if interrupted else 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
