import numpy as np

from chebyshev_hofbauer_resonances.general_tent_map.adjacency_to_super import (
    create_partial_super_adjacency,
)
from chebyshev_hofbauer_resonances.general_tent_map.hofbauer_tower import (
    create_adjacency_matricies,
)

from chebyshev_hofbauer_resonances.general_tent_map.ulams_method import ulams_method


def construct_transfer_operators(inverses, derivatives):
    """
    Constructs the transfer operator for each "segemnt" of the piecewise function.

    Parameters
    ----------
    inverses : list
        The list of inverse functions for each segment.
    derivatives : list
        The list of derivative functions for each segment.

    Returns
    -------
    transfer_operators : list
        The list of transfer operator functions for each segment.
    """

    return [
        lambda phi, i=i: lambda x: phi(inverses[i](x))
        / abs(derivatives[i](inverses[i](x)))
        for i in range(len(inverses))
    ]


def create_super_adjacency(domains, adj_matrices, transfer_operators, N, K, depth):
    """
    Create the super adjacency matrix from the adjacency matrices and transfer operators.

    Parameters
    ----------
    domains : list
        The list of domains.
    adj_matrices : list
        The list of adjacency matrices.
    transfer_operators : list
        The list of transfer operator functions.
    K : integer
        The order of the Chebyshev nodes, taken to be the order of the DCT.
    N : integer
        The order of the Chebyshev polynomials to use.
    depth : int
        The depth of the approximation.
    Returns
    -------
    super_adjacency : ndarray
        The super adjacency matrix.
    """
    partial_super_adjacency = [
        create_partial_super_adjacency(
            adj_matrices[i], domains, transfer_operators[i], K, N, depth
        )
        for i in range(len(adj_matrices))
    ]

    super_adjacency = np.sum(partial_super_adjacency, axis=0)

    return super_adjacency


def approx_super_adjacency(
    function_domains, functions, inverses, derivatives, N, K, depth
):
    """
    Create the super adjacency matrix approximation for the given piecewise function.
    Parameters
    ----------
    function_domains : list
        The list of domains for each segment of the piecewise function.
    functions : list
        The list of functions for each segment of the piecewise function.
    inverses : list
        The list of inverse functions for each segment.
    derivatives : list
        The list of derivative functions for each segment.
    K : integer
        The order of the Chebyshev nodes, taken to be the order of the DCT.
    N : integer
        The order of the Chebyshev polynomials to use.
    depth : int
        The depth of the approximation.
    Returns
    -------
    super_adjacency : ndarray
        The super adjacency matrix approximation.
    """
    domains, adj_matrices = create_adjacency_matricies(
        function_domains, functions, depth=depth
    )
    transfer_operators = construct_transfer_operators(inverses, derivatives)
    super_adjacency = create_super_adjacency(
        domains, adj_matrices, transfer_operators, N, K, 1
    )

    return super_adjacency


def approx_ulams(function_domains, functions, inverses, derivatives, N, M):
    """
    Create the Ulam's method approximation for the given piecewise function.
    Parameters
    ----------
    function_domains : list
        The list of domains for each segment of the piecewise function.
    functions : list
        The list of functions for each segment of the piecewise function.
    inverses : list
        The list of inverse functions for each segment.
    derivatives : list
        The list of derivative functions for each segment.
    N : integer
        The number of Ulam bins.
    M : integer
        The number of samples per Ulam bin.
    Returns
    -------
    L : ndarray
        The Ulam's method approximation of the transfer operator.
    """

    def f(x):
        for i, domain in enumerate(function_domains):
            if domain[0] <= x <= domain[1]:
                return functions[i](x)
        raise ValueError(f"x={x} is not in any domain")

    f = np.vectorize(f)

    return ulams_method(N, M, f)
