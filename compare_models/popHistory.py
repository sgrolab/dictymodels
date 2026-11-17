"""
What a population of cells keeps of a run, how it is given its noise, and what the compiled
loops of the five models share.

A run steps N cells through len(t) time points. Every variable of every cell at every time
point is N x len(t) numbers, which for a sweep job is gigabytes, so by default run() keeps
only the mean of each variable over the cells, as <name>_mean. Call run(..., keep_cells=True)
to also keep each variable for every cell, as <name>, an N x len(t) array. That is what the
models always did before and is meant for looking inside a run. Variables that all cells
share, such as the external cAMP of the Maeda, Gregor and Sgro models, are always kept.

The noise r is either the whole len(t) x N array, as before, or the numpy random generator
to draw it from. Given a generator, a run draws the noise of each time step when it gets
there, the numbers the whole array would have held, and keeps none of it.

The loop over the time steps is the method steps() of each class, which numba compiles the
first time it is called; that makes it about ten times faster. Its statements are the numpy
statements the loop always had.

numba keeps each compiled loop in compare_models/__pycache__ and uses it again in later
runs, which saves most of a second each time, until the model file that holds the loop is
changed. It does not notice a change in this file: after changing a compiled function
here, delete that folder.

A compiled loop gives the same numbers as the same statements run by Python, to the last
digit, if three things are kept to when changing one:
  - use np_mean and np_sum, below, for np.mean and np.sum. numba's own add the elements up
    one after another, while numpy adds in eight interleaved partial sums, and the two
    round differently;
  - write a power other than a square with a decimal point, A**3.0, not A**3. numpy takes
    every power but a square by the pow of the C library; numba does so only for an
    exponent that is not an integer, and multiplies otherwise;
  - leave fastmath off, and keep error_model='numpy', with which a division by zero gives
    inf or nan as in numpy instead of stopping the run.
numba compiles for the processor it runs on; other machines can differ in the last digit.
"""
import math

import numpy as np
from numba import njit


def noise_source(r, steps, N, dt, sigma=1):
    """Return the noise of a run of N cells as a compiled loop takes it: see noise().

    r is the noise itself, an array of time steps x cells, or a numpy random generator.
    From a generator the loop draws the noise of each time step as
    math.sqrt(dt)*r.normal(0,1,N)*sigma. Those are the numbers that
    math.sqrt(dt)*r.normal(0,1,size=(steps,N))*sigma would hold, which is how the scripts
    used to make the array; sigma is used only then."""
    if isinstance(r, np.random.Generator):
        r.normal(0, 1, N)   # the row of the first time step, which takes no noise
        return np.zeros((0, N)), r, math.sqrt(dt), float(sigma)
    if callable(r):
        raise TypeError('the noise must be an array of time steps x cells or a numpy random '
                        'generator, not a function')
    rows = np.asarray(r, dtype=float)
    if rows.ndim == 0 or len(rows) < steps:
        raise ValueError(f'the noise needs a row for each of the {steps} time steps, but '
                         f'its shape is {rows.shape}')
    rows = rows.reshape(len(rows), -1)
    if rows.shape[1] != N:      # one value per time step is the same for every cell
        rows = np.broadcast_to(rows, (len(rows), N))
    return np.ascontiguousarray(rows), np.random.default_rng(), 0.0, 0.0


@njit(error_model='numpy')
def noise(r, i):
    """Return the noise of the cells at time step i, r being what noise_source() made of
    the noise given to run(): row i of the array, or drawn now from the generator."""
    rows, generator, sqrt_dt, sigma = r
    if len(rows) > 0:
        return rows[i]
    return sqrt_dt*generator.normal(0.0, 1.0, rows.shape[1])*sigma


def stimulus(values, steps):
    """Return the stimulus, with a value for each time step, as a compiled loop takes it."""
    values = np.ascontiguousarray(values, dtype=float)
    if values.ndim == 0 or len(values) < steps:
        raise ValueError(f'the stimulus needs a value for each of the {steps} time steps, '
                         f'but its shape is {values.shape}')
    return values


@njit(error_model='numpy')
def block_sum(a, start, n):
    """Return the sum of a[start:start + n], n at most 128, in the order numpy adds: in
    eight interleaved partial sums, then what is left over."""
    if n < 8:
        res = -0.0
        for i in range(n):
            res += a[start + i]
        return res
    r0, r1, r2, r3 = a[start], a[start + 1], a[start + 2], a[start + 3]
    r4, r5, r6, r7 = a[start + 4], a[start + 5], a[start + 6], a[start + 7]
    i = 8
    while i < n - n % 8:
        r0 += a[start + i]
        r1 += a[start + i + 1]
        r2 += a[start + i + 2]
        r3 += a[start + i + 3]
        r4 += a[start + i + 4]
        r5 += a[start + i + 5]
        r6 += a[start + i + 6]
        r7 += a[start + i + 7]
        i += 8
    res = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
    while i < n:
        res += a[start + i]
        i += 1
    return res


@njit(error_model='numpy')
def pairwise_sum(a, start, n):
    """Return the sum of a[start:start + n], added up in the order in which numpy adds.

    numpy halves a vector of more than 128, and each half again, with a function that
    calls itself. Here the halves still to be added wait on a stack instead, because numba
    cannot keep compiled code that calls itself for a later run: loading it crashes."""
    if n <= 128:
        return block_sum(a, start, n)
    starts = np.empty(64, np.int64)
    sizes = np.empty(64, np.int64)
    stages = np.empty(64, np.int64)
    lefts = np.empty(64)
    top = 0
    starts[0], sizes[0], stages[0] = start, n, 0
    result = 0.0
    while top >= 0:
        s, m = starts[top], sizes[top]
        if m <= 128:
            result = block_sum(a, s, m)
            top -= 1
            continue
        half = m // 2
        half -= half % 8
        if stages[top] == 0:        # first the left half
            stages[top] = 1
            top += 1
            starts[top], sizes[top], stages[top] = s, half, 0
        elif stages[top] == 1:      # then the right half
            lefts[top] = result
            stages[top] = 2
            top += 1
            starts[top], sizes[top], stages[top] = s + half, m - half, 0
        else:                       # then the two together
            result = lefts[top] + result
            top -= 1
    return result


@njit(error_model='numpy')
def np_sum(a):
    """np.sum(a) of a vector a for a compiled loop, the same to the last digit."""
    total = 0.0
    for start in range(0, len(a), 8192):    # numpy adds a long vector up 8192 at a time
        total += pairwise_sum(a, start, min(8192, len(a) - start))
    return total


@njit(error_model='numpy')
def np_mean(a):
    """np.mean(a) of a vector a for a compiled loop, the same to the last digit."""
    return np_sum(a)/len(a)


@njit(error_model='numpy')
def np_heaviside(x1, x2):
    """np.heaviside(x1, x2) of a vector x1 and a number x2, which numba does not have."""
    out = np.empty(len(x1))
    for k in range(len(x1)):
        if x1[k] < 0:
            out[k] = 0.0
        elif x1[k] == 0:
            out[k] = x2
        elif x1[k] > 0:
            out[k] = 1.0
        else:               # not a number
            out[k] = x1[k]
    return out


@njit(error_model='numpy')
def put(history, i, *values):
    """Keep the variables of time step i in history, the arrays of a History. values are
    the variables in the order of its names, each a vector with one value for each cell.

    The mean adds the cells up one after another, because that is how np.mean(x, axis=0)
    treats an x of cells x time steps, which the scripts used to take; np_mean rounds
    differently."""
    means, cells = history
    for k in range(len(values)):
        value = values[k]
        total = -0.0
        for cell in range(len(value)):
            total += value[cell]
        means[k, i] = total/len(value)
        if cells.shape[1] > 0:
            for cell in range(len(value)):
                cells[k, cell, i] = value[cell]


class History:
    """What is kept of the variables of the cells of a population during its run()."""

    def __init__(self, pop, names, first, keep_cells):
        """pop is the population, names the variables in the order put() is given them, and
        first their values at the first time step, one for each cell."""
        self.pop = pop
        self.names = names
        self.keep_cells = keep_cells
        steps = len(pop.t)
        means = np.zeros([len(names), steps])
        cells = np.zeros([len(names), len(first[0]), steps] if keep_cells else [0, 0, 0])
        self.arrays = (means, cells)    # what put() fills, also from inside a compiled loop
        put(self.arrays, 0, *first)

    def finish(self):
        """Give the population what was kept: <name>_mean, and <name> if cells were kept."""
        means, cells = self.arrays
        for k, name in enumerate(self.names):
            setattr(self.pop, name + '_mean', means[k])
            if self.keep_cells:
                setattr(self.pop, name, cells[k])


class CellPopulation:
    """Base of the population classes whose cells differ; explains a variable that is absent."""

    def __getattr__(self, name):
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'. The variables of the "
            "cells exist only after run(): their means over the cells as <name>_mean, and "
            "with run(..., keep_cells=True) also the values of every cell as <name>.")
