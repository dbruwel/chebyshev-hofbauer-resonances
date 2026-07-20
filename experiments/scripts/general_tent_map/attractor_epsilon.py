# imports
import concurrent.futures

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from tqdm import tqdm

from chebyshev_hofbauer_resonances.general_tent_map.approx_transfer_op import (
    approx_super_adjacency,
)


# function definition
def a_map(eps):
    function_domains = [(0, 0.5), (0.5, 1)]
    functions = [lambda x: eps + x * (1 + (2 - 4 * eps) * x), lambda x: 2 * x - 1]
    inverses = [
        lambda y: (
            (-1 + np.sqrt((1 - 4 * eps) ** 2 + (8 - 16 * eps) * y)) / (4 - 8 * eps)
        ),
        lambda y: (1 + y) / 2,
    ]
    derivatives = [lambda x: 1 + (4 - 8 * eps) * x, lambda _: 2]

    return function_domains, functions, inverses, derivatives


def get_lam2(eps):
    function_domains, functions, inverses, derivatives = a_map(eps)

    try:
        super_adjacency = approx_super_adjacency(
            function_domains,
            functions,
            inverses,
            derivatives,
            N=50,
            K=50,
            depth=100,
        )

        sparse_adj = sp.csr_matrix(super_adjacency)
        evals_super_adj = spla.eigs(
            sparse_adj, k=2, which="LM", return_eigenvectors=False
        )
        evals_super_adj = evals_super_adj[np.argsort(-np.abs(evals_super_adj))]
        lam2 = evals_super_adj[1]
        return np.abs(lam2)

    except Exception:
        print(f"Failed for epsilon: {eps}")
        return 0


if __name__ == "__main__":
    lam2s = []
    epss = np.linspace(0, 0.75, 750_001)

    NUM_CPUS = 50

    ress = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_CPUS) as executor:
        results_generator = executor.map(get_lam2, epss)
        for res in tqdm(results_generator, total=len(epss), desc="Processing Epsilons"):
            ress.append(res)

    df = pd.DataFrame(ress, index=epss, columns=["lam_2"])
    df.to_csv("res.csv")
