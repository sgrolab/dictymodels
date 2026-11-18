"""
What a population of cells keeps of a run, and how it is given its noise.

A run steps N cells through len(t) time points. Every variable of every cell at every time
point is N x len(t) numbers, which for a sweep job is gigabytes, so by default run() keeps
only the mean of each variable over the cells, as <name>_mean. Call run(..., keep_cells=True)
to also keep each variable for every cell, as <name>, an N x len(t) array. That is what the
models always did before and is meant for looking inside a run. Variables that all cells
share, such as the external cAMP of the Maeda, Gregor and Sgro models, are always kept.

The noise r is either the whole len(t) x N array, as before, or a function that returns the
next rows of it when called with a number of rows. Given a function, a run draws its noise as
it goes and never holds more than BLOCK rows.
"""
import numpy as np

BLOCK = 1000    # time steps of noise, and of the variables not yet averaged, held at once


def noise_rows(r, steps):
    """Iterate over r[0], r[1], ..., the noise of each cell at each of the time steps.

    r is a steps x N array, or a function r(n) that returns the next n rows of one."""
    if not callable(r):
        yield from r
        return
    for start in range(0, steps, BLOCK):
        yield from r(min(BLOCK, steps - start))


def mean_over_cells(rows):
    """Return the mean of each row of rows, which is time steps x cells.

    The cells are added up in order because that is how np.mean(x, axis=0) treats an x of
    cells x time steps, so the result is the same to the last digit. np.mean(rows, axis=1)
    adds them in pairs and rounds differently."""
    total = rows[:, 0].copy()
    for cell in range(1, rows.shape[1]):
        total += rows[:, cell]
    return total / rows.shape[1]


class History:
    """Collects the variables of the cells of a population as its run() steps through time."""

    def __init__(self, pop, names, first, keep_cells):
        """pop is the population, names the variables in the order put() is given them, and
        first their values at the first time step, one for each cell."""
        self.pop = pop
        self.names = names
        self.keep_cells = keep_cells
        self.steps = len(pop.t)
        N = len(first[0])
        if keep_cells:
            self.cells = [np.zeros([N, self.steps]) for name in names]
        else:
            self.means = [np.zeros(self.steps) for name in names]
            self.blocks = [np.zeros([BLOCK, N]) for name in names]
        self.put(0, *first)

    def put(self, i, *values):
        """Take the variables at time step i. The steps must come in order."""
        if self.keep_cells:
            for cells, value in zip(self.cells, values):
                cells[:, i] = value
            return
        row = i % BLOCK
        for block, value in zip(self.blocks, values):
            block[row] = value
        if row == BLOCK - 1 or i == self.steps - 1:
            for mean, block in zip(self.means, self.blocks):
                mean[i - row:i + 1] = mean_over_cells(block[:row + 1])

    def finish(self):
        """Give the population what was kept: <name>_mean, and <name> if cells were kept."""
        for k, name in enumerate(self.names):
            if self.keep_cells:
                setattr(self.pop, name, self.cells[k])
                setattr(self.pop, name + '_mean', mean_over_cells(self.cells[k].T))
            else:
                setattr(self.pop, name + '_mean', self.means[k])


class CellPopulation:
    """Base of the population classes whose cells differ; explains a variable that is absent."""

    def __getattr__(self, name):
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'. The variables of the "
            "cells exist only after run(): their means over the cells as <name>_mean, and "
            "with run(..., keep_cells=True) also the values of every cell as <name>.")
