# imports
import concurrent.futures

import numpy as np
import pandas as pd
from tqdm import tqdm

from chebyshev_hofbauer_resonances.general_tent_map.approx_transfer_op import (
    approx_super_adjacency,
)


# function definition
def t_map(alpha):
    function_domains = [(0, 0.5), (0.5, 1)]
    functions = [lambda x: alpha * x, lambda x: alpha * (1 - x)]
    inverses = [lambda y: y / alpha, lambda y: 1 - (y / alpha)]
    derivatives = [lambda _: alpha, lambda _: -alpha]

    return function_domains, functions, inverses, derivatives


def get_lam2(alpha):
    function_domains, functions, inverses, derivatives = t_map(alpha)

    super_adjacency = approx_super_adjacency(
        function_domains,
        functions,
        inverses,
        derivatives,
        N=2,
        K=2,
        depth=50,
    )

    evals_super_adj = np.linalg.eigvals(super_adjacency)
    evals_super_adj = evals_super_adj[np.argsort(-np.abs(evals_super_adj))]

    lam2 = evals_super_adj[1]
    return lam2


if __name__ == "__main__":
    lam2s = []
    alphas = np.linspace(1.4, 2, 600_001)

    NUM_CPUS = 6

    ress = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_CPUS) as executor:
        results_generator = executor.map(get_lam2, alphas)
        for res in tqdm(results_generator, total=len(alphas), desc="Processing Alphas"):
            ress.append(res)

    df = pd.DataFrame(ress, index=alphas, columns=["lam_2"])
    df.to_csv("res.csv")
