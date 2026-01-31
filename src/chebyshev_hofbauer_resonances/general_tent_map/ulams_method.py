import matplotlib.pyplot as plt
import numpy as np

"""
Old code repeated.
"""


def ulams_method(N, M, f):
    """
    Compute Ulam's method approximation of the transfer operator.

    Parameters
    ----------
    N : int
        The number of bins to partition the interval [0, 1].
    M : int
        The number of sample points per bin.
    f : callable
        The map function to approximate the transfer operator for.

    Returns
    -------
    L : ndarray
        The Ulam's method approximation matrix of shape (N, N).
    """
    bins = np.linspace(0, 1, N + 1)

    L = np.zeros((N, N))

    for i in range(N):
        x_samples = np.linspace(bins[i], bins[i + 1], M)

        x_next = f(x_samples)

        bin_indices = np.digitize(x_next, bins) - 1
        bin_indices = np.clip(bin_indices, 0, N - 1)

        for j in bin_indices:
            L[i, j] += 1

    L /= L.sum(axis=1, keepdims=True)

    return L
