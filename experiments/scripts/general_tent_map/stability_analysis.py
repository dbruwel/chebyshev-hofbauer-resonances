# imports
import numpy as np
import pyperf

from chebyshev_hofbauer_resonances.general_tent_map.approx_transfer_op import (
    approx_super_adjacency,
    approx_ulams,
)

# tent map definition
alpha = 1.2
function_domains = [(0, 0.5), (0.5, 1)]
functions = [lambda x: 1.2 * x, lambda x: 1.2 * (1 - x)]
inverses = [lambda y: y / 1.2, lambda y: 1 - (y / 1.2)]
derivatives = [lambda _: 1.2, lambda _: -1.2]


# functions
def super_adj_error(N, K, depth):
    super_adjacency = approx_super_adjacency(
        function_domains,
        functions,
        inverses,
        derivatives,
        N=N,
        K=K,
        depth=depth,
    )

    evals_super_adj = np.linalg.eigvals(super_adjacency)
    evals_super_adj = evals_super_adj[np.argsort(-np.abs(evals_super_adj))]

    max_eval = np.abs(evals_super_adj).max()

    error = np.abs(1 - max_eval)

    return error


def ulam_error(N, M):
    L_ulam = approx_ulams(
        function_domains,
        functions,
        inverses,
        derivatives,
        N=N,
        M=M,
    )

    evals_ulam = np.linalg.eigvals(L_ulam)
    evals_ulam = evals_ulam[np.argsort(-np.abs(evals_ulam))]

    max_eval = np.abs(evals_ulam).max()

    error = np.abs(1 - max_eval)

    return error


if __name__ == "__main__":
    runner = pyperf.Runner()

    N = 10
    K = 10
    depth = 5
    runner.bench_func("super_adj_error", lambda: super_adj_error(N, K, depth))
