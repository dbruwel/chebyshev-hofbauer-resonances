'''import numpy as np
from scripts.operator_approx import cheb_op_ap


def generate_opp_approx(i, j, domains, generate, L, N, K):
    """
    Takes the i, j element of the adjacency matrix and generates the operator approximation if generate is True.
    Otherwise returns the 0 matrix. 

    Parameters
    ----------
    i : int
        The row index of the adjacency matrix.
    j : int
        The column index of the adjacency matrix.
    domains : list
        The list of domains.
    generate : bool
        Whether to generate the operator approximation.
    L : function
        The operator to approximate.
    K : integer
        The order of the Chebyshev nodes, taken to be the order of the DCT.
    N : integer
        The order of the Chebyshev polynomials to use.

    Returns
    -------
    L_hat : ndarray
        The operator approximation.
    """
    if generate:
        # There is a good chance this is backwards because I can't remember the order of adjacency matrices.
        final_domain = domains[i]
        initial_domain  = domains[j]

        L_hat_T = cheb_op_ap(L, K, N, initial_domain=initial_domain, final_domain=final_domain)
        L_hat = L_hat_T.T
        return L_hat
    else:
        return np.zeros((K, N))
    
def create_super_adjacency(adjacency_matrix, domains, L, N, K):
    """
    Create the super adjacency matrix from the adjacency matrix.
    
    Parameters
    ----------
    adjacency_matrix : ndarray
        The adjacency matrix.
    domains : list
        The list of domains.
    L : function
        The operator to approximate.
    K : integer
        The order of the Chebyshev nodes, taken to be the order of the DCT.
    N : integer
        The order of the Chebyshev polynomials to use.

    Returns
    -------
    super_adjacency : ndarray
        The super adjacency matrix.
    """
    super_adjacency = np.empty((adjacency_matrix.shape[0], adjacency_matrix.shape[1]), dtype=object)

    for i in range(adjacency_matrix.shape[0]):
        for j in range(adjacency_matrix.shape[1]):
            super_adjacency[i, j] = generate_opp_approx(i, j, domains, adjacency_matrix[i, j], L, N, K)

    super_adjacency = np.block(super_adjacency.tolist())

    return super_adjacency'''



import numpy as np
from scripts.operator_approx import cheb_op_ap

def generate_opp_approx(i, j, domains, generate, L, N, K, depth):
    """
    Takes the i, j element of the adjacency matrix and generates the operator approximation if generate is True.
    Otherwise returns the 0 matrix. 

    Parameters
    ----------
    i : int
        The row index of the adjacency matrix.
    j : int
        The column index of the adjacency matrix.
    domains : list
        The list of domains.
    generate : bool
        Whether to generate the operator approximation.
    L : function
        The operator to approximate.
    K : integer
        The order of the Chebyshev nodes, taken to be the order of the DCT.
    N : integer
        The order of the Chebyshev polynomials to use.
    depth : int
        The depth of the approximation.

    Returns
    -------
    L_hat : ndarray
        The operator approximation.
    """
    if generate:
        # There is a good chance this is backwards because I can't remember the order of adjacency matrices.
        final_domain = domains[i]
        initial_domain = domains[j]

        # Pass depth to cheb_op_ap
        L_hat_T = cheb_op_ap(L, K, N, initial_domain=initial_domain, final_domain=final_domain, depth=depth)
        L_hat = L_hat_T.T
        return L_hat
    else:
        return np.zeros((K, N))
    
def create_super_adjacency(adjacency_matrix, domains, L, N, K, depth):
    """
    Create the super adjacency matrix from the adjacency matrix.
    
    Parameters
    ----------
    adjacency_matrix : ndarray
        The adjacency matrix.
    domains : list
        The list of domains.
    L : function
        The operator to approximate.
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
    super_adjacency = np.empty((adjacency_matrix.shape[0], adjacency_matrix.shape[1]), dtype=object)

    for i in range(adjacency_matrix.shape[0]):
        for j in range(adjacency_matrix.shape[1]):
            # Pass depth to generate_opp_approx
            super_adjacency[i, j] = generate_opp_approx(i, j, domains, adjacency_matrix[i, j], L, N, K, depth)

    super_adjacency = np.block(super_adjacency.tolist())

    return super_adjacency
