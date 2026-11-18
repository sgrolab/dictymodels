#!/usr/bin/env python
"""Run the simulation scripts that generate the data for the figures.

Each script is run in its own Python process, exactly as it would be by hand. As many run
at a time as the computer has processors, or as --jobs says. The scripts of the density and
flow-rate sweep simulate one grid point per run, so each of them is run once for every point
of its grid; that is 9,150 runs and nearly all of the work. Starting Python and loading its
libraries would be a large part of so short a run, so the sweep is run by processes that
stay: each runs a sweep script for one grid point, from its first line to its last, and
then takes the next grid point that is left. When every grid point of the sweep has its
file, the consolidate script collects them for the figures.

Without arguments every script is run except those in EXCLUDE; name scripts, or the start
of the names of several, to run a subset. The data directory is the one the scripts use
(ROOT_DIR from the environment or .env).

Ctrl-C stops a run, as does SIGTERM, which is what `kill` sends; --resume then carries on
with the grid points that have no file yet.
"""
import argparse
import contextlib
import gc
import io
import json
import os
import queue
import runpy
import signal
import subprocess
import sys
import threading
import time
import traceback
import warnings
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
SWEEP_PREFIX = 'Population_DensityFlowSweep_'
CONSOLIDATE = SWEEP_PREFIX + 'Consolidate.py'

# Scripts left out of a full run, and why; they can still be named on the command line
EXCLUDE = {'Population_add_cAMP.py': 'no figure uses its results'}

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
    found = set(named(SWEEP_PREFIX + '*.py')) - {CONSOLIDATE}
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
        return len(os.sched_getaffinity(0))     # Linux: no more than taskset allows
    except AttributeError:                      # macOS and Windows cannot say
        return os.cpu_count() or 1


def hms(seconds):
    seconds = int(round(seconds))
    return f'{seconds // 3600}:{seconds // 60 % 60:02d}:{seconds % 60:02d}'


def run_points():
    """What `runSimulations.py --points` does, for Runner.run_point(): run sweep scripts for
    one grid point after another, here in this process. Each line of its input asks for a
    run, as SCRIPT I J, and it ends when there are no more.

    Each run is the whole script, from its first line, as `python SCRIPT I J` runs it. Only
    the libraries, and the loops that numba has compiled, stay from one run to the next.
    Says how each run went in a line of JSON."""
    report = sys.stdout
    try:
        for line in sys.stdin:
            if not line.strip():
                continue
            script, *indices = line.split()
            path = str(CODE_DIR / script)
            text, code = io.StringIO(), 0
            own_argv, sys.argv = sys.argv, [path] + indices
            try:
                with contextlib.redirect_stdout(text), contextlib.redirect_stderr(text), \
                        warnings.catch_warnings():  # so that a warning is given in every run
                    runpy.run_path(path, run_name='__main__')
            except SystemExit as stop:
                if stop.code is not None and not isinstance(stop.code, int):
                    text.write(f'{stop.code}\n')
                code = stop.code if isinstance(stop.code, int) else int(stop.code is not None)
            except Exception:
                code = 1
                text.write(traceback.format_exc())
            finally:
                sys.argv = own_argv
            gc.collect()    # the arrays of the run that has ended
            print(json.dumps({'done': [script] + indices, 'code': code,
                              'text': text.getvalue()}), file=report, flush=True)
    except KeyboardInterrupt:
        return 130
    return 0


class Runner:
    """Makes the runs in the queue work, in processes, and stops them when told to. How each
    run went is put in the queue results as soon as it is known."""

    def __init__(self, log_dir, stop_on_error):
        self.log_dir = log_dir
        self.stop_on_error = stop_on_error
        self.stop = threading.Event()   # set to start no more runs
        self.lock = threading.Lock()
        self.processes = set()
        self.work = queue.SimpleQueue()     # of each run to make: script, indices, file
        self.results = queue.Queue()

    def start(self, command, **options):
        """Start a process unless no more are to be started; returns it, or None."""
        with self.lock:
            if self.stop.is_set():
                return None
            process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True,
                                       stderr=subprocess.STDOUT, errors='replace', **options)
            self.processes.add(process)
        return process

    def end(self, process):
        """Wait for a process that is ending; returns its exit status."""
        process.wait()
        with self.lock:
            self.processes.discard(process)
        return process.returncode

    def finish(self, result, output=None):
        """Take in how a run went: 'code' is its exit status. output is the file that the
        run of a grid point must have written."""
        script, indices = result['script'], result['indices']
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
        self.results.put(result)
        return result

    def run(self, script):
        """Run a script in a process of its own; returns how it went, or None if it was not
        started."""
        start = time.monotonic()
        process = self.start([sys.executable, str(CODE_DIR / script)])
        if process is None:
            return None
        try:
            text = process.communicate()[0]
        except KeyboardInterrupt:   # in the main thread, which runs the consolidate script
            process.terminate()
            self.end(process)       # it has ended before this is passed on
            raise
        return self.finish({'script': script, 'indices': (), 'code': self.end(process),
                            'seconds': time.monotonic() - start, 'text': text})

    def run_point(self, process, script, indices, output):
        """Run a sweep script for one grid point, whose file is output, in a process that
        stays for the next: see run_points(). The process is started if None is given.
        Returns it, or None if it died, in which case this grid point has failed."""
        start = time.monotonic()
        if process is None:
            process = self.start([sys.executable, str(CODE_DIR / 'runSimulations.py'),
                                  '--points'], stdin=subprocess.PIPE)
            if process is None:
                return None
        try:
            process.stdin.write(' '.join([script] + [str(i) for i in indices]) + '\n')
            process.stdin.flush()
        except OSError:     # it has died: there will be nothing to read either
            pass
        text, said = '', None
        for line in iter(process.stdout.readline, ''):
            at = line.find('{"done"')   # not always where the line starts: what a run wrote
            if at >= 0:                 # past Python may not have ended its line
                try:
                    said = json.loads(line[at:])
                    text += (line[:at] + '\n') if at else ''
                    break
                except ValueError:
                    pass
            text += line    # what the run wrote past Python, or the process itself
        if said is None:    # the process ended before the run did
            code, process = self.end(process) or 'died', None
        else:
            code, text = said['code'], text + said['text']
        self.finish({'script': script, 'indices': indices, 'code': code,
                     'seconds': time.monotonic() - start, 'text': text}, output)
        return process

    def take(self):
        """What each thread does: make the runs in the queue, one after another, until none
        is left. Its grid points all go to one process, which takes the place of the
        Python that each would otherwise start."""
        process = None
        try:
            while not self.stop.is_set():
                try:
                    script, indices, output = self.work.get_nowait()
                except queue.Empty:
                    break
                if indices:
                    process = self.run_point(process, script, indices, output)
                else:
                    self.run(script)
        finally:
            if process is not None:
                process.stdin.close()   # with no more to read, it ends
                self.end(process)

    def interrupt(self):
        """Start no more runs and end those that are under way."""
        self.stop.set()
        with self.lock:
            for process in self.processes:
                process.terminate()


def stop_as_ctrl_c(signum, frame):
    """What SIGTERM does to a run: the same as Ctrl-C, so that the runs under way are
    stopped and their unfinished files removed."""
    raise KeyboardInterrupt


def show_output(text, lines=20):
    """Print the end of what a run wrote, indented."""
    for line in text.rstrip().split('\n')[-lines:] if text.strip() else []:
        print('    ' + line)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv == ['--points']:    # one of the processes of Runner.run_point()
        return run_points()
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
    jobs = max(1, min(args.jobs, len(once) + total))    # no more than there is to run

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
    signal.signal(signal.SIGTERM, stop_as_ctrl_c)
    began = time.monotonic()
    reported = began
    interrupted = False
    print()
    results = []
    for s in once:
        runner.work.put((s, (), None))
    for s, points in todo.items():
        for indices, output in points:
            runner.work.put((s, indices, output))
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(runner.take) for _ in range(jobs)]
        for future in futures:      # the queue also says when a thread has ended
            future.add_done_callback(runner.results.put)
        try:
            done = failed = 0
            while futures:
                try:
                    result = runner.results.get(timeout=1)
                except queue.Empty:
                    continue
                if isinstance(result, Future):
                    futures.remove(result)
                    result.result()     # pass on what went wrong in the thread, if anything
                    continue
                results.append(result)
                name = ' '.join([result['script']] + [str(i) for i in result['indices']])
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
    cut = set()     # the scripts with a run that the interrupt stopped
    while not runner.results.empty():   # all threads have ended: what an interrupt left
        result = runner.results.get()
        if isinstance(result, Future):
            continue
        if interrupted and result['code'] != 0:
            cut.add(result['script'])   # stopped, not failed: its file is gone, to run again
        else:
            results.append(result)

    # the consolidate script, if it was asked for and the sweep is all there
    consolidated = None
    started = False
    if CONSOLIDATE in scripts and not runner.stop.is_set():
        try:
            missing = {
                s: sum(not output.is_file() for _, output in grid_points(s, PROJECT_DIR))
                for s in all_scripts() if s.startswith(SWEEP_PREFIX) and s != CONSOLIDATE}
            if any(missing.values()):
                consolidated = 'not run: no file for ' + ', '.join(
                    f'{n:,} grid point(s) of {s}' for s, n in missing.items() if n)
            else:
                started = True
                consolidated = runner.run(CONSOLIDATE)
        except KeyboardInterrupt:
            interrupted = True
            consolidated = 'interrupted before it ' + \
                           ('had ended' if started else 'was started')
            print('\ninterrupted: stopping the consolidate script' if started else
                  '\ninterrupted: the consolidate script is not started')
            runner.interrupt()

    print("\n===== summary (the time of a sweep script is that of its runs added up) =====")
    bad = 0
    for s in scripts:
        ran = [r for r in results if r['script'] == s]
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
                state = 'stopped' if started else 'not run'
                print(f'{state}          {s}')
                if consolidated:
                    print(f'    {consolidated}')
        elif s in sweep:
            ok = not failures and len(ran) == len(todo[s])
            line = f'{len(ran) - len(failures):,} of {len(todo[s]):,} runs ok'
            if failures:
                failures.sort(key=lambda r: r['indices'])
                line += '; failed: ' + ', '.join(
                    ' '.join(str(i) for i in r['indices']) for r in failures[:10])
                line += ', ...' if len(failures) > 10 else ''
            state = 'FAILED' if failures else 'stopped' if ran or s in cut else 'not run'
            print(f'{"ok    " if ok else state}  {seconds}  {s}  {line}')
        else:
            ok = bool(ran) and not failures
            state = 'FAILED' if ran else 'stopped' if s in cut else 'not run'
            print(f'{"ok    " if ok else state}  {seconds}  {s}')
        bad += not ok
    print(f'{len(scripts) - bad} script(s) ok, {bad} not; '
          f'{hms(time.monotonic() - began)} elapsed')
    return 130 if interrupted else 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
