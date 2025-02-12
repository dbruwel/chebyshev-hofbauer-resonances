import numpy as np
import matplotlib.pyplot as plt

def ulams_method(N, M, f):
    bins = np.linspace(0, 1, N+1)

    L = np.zeros((N, N))

    for i in range(N):
        x_samples = np.random.uniform(bins[i], bins[i+1], M)
        
        x_next = f(x_samples)
        
        bin_indices = np.digitize(x_next, bins) - 1
        bin_indices = np.clip(bin_indices, 0, N-1)
        
        for j in bin_indices:
            L[i, j] += 1

    L /= L.sum(axis=1, keepdims=True)

    return L