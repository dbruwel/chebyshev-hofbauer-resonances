import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.chebyshev import chebfit
from scipy.fftpack import dct
from scipy.special import chebyt


def cheb_op_ap(L, K, N, initial_domain=(-1, 1), final_domain=(-1, 1), depth=1):
    """
    Return the Chebyshev matrix approximation of an operator with a depth refinement.

    Parameters
    ----------
    L : function
        The input operator. Should be a function that takes in a function and returns a function.
    K : integer
        The order of the Chebyshev nodes, taken to be the order of the DCT.
    N : integer
        The order of the Chebyshev polynomials to use.
    initial_domain : tuple, optional
        The initial domain of the operator. The default is (-1, 1).
    final_domain : tuple, optional
        The final domain of the operator. The default is (-1, 1).
    depth : integer, optional
        The depth of refinement for the approximation. Default is 1.

    Returns
    -------
    L_hat : ndarray
        The matrix approximation of the operator, refined based on depth.
    """
    assert N == K, "currently only works for K = N"

    k = np.arange(0, K)
    theta = np.pi * (2 * k + 1) / (2 * K)
    x = np.cos(theta)
    x = linear_map(x, final_domain)

    y = np.array([L(domain_restricted_chebyt(n, initial_domain))(x) for n in range(N)])

    L_hat = dct(y, type=2, axis=1) / y.shape[1]
    L_hat[:, 0] = L_hat[:, 0] / 2

    # Apply depth refinement (for example, running refinement iterations here)
    if depth > 1:
        for _ in range(depth - 1):
            # Recalculate approximation by applying some refinement process (example: further smoothing)
            L_hat = dct(L_hat, type=2, axis=1) / L_hat.shape[1]
            L_hat[:, 0] = L_hat[:, 0] / 2

    return L_hat


def linear_map(values, domain):
    """
    Linearly map values from the domain [-1, 1] to the domain [a, b].

    Parameters
    ----------
    values : ndarray
        The input array of values in the domain [-1, 1].
    domain : tuple
        The target domain to map the

    Returns
    -------
    mapped_values : ndarray
        The array of values mapped to the domain [a, b].

    Examples
    --------
    >>> values = np.array([-1, -0.5, 0, 0.5, 1])
    >>> a, b = 2, 3
    >>> linear_map(values, (a, b))
    array([2. , 2.25, 2.5 , 2.75, 3. ])
    """
    a, b = domain
    return a + (values + 1) * (b - a) / 2


def inverse_linear_map(values, domain):
    """
    Linearly map values from the domain [a, b] to the domain [-1, 1].

    Parameters
    ----------
    values : ndarray
        The input array of values in the domain [a, b].
    domain : tuple
        The target domain to map the

    Returns
    -------
    mapped_values : ndarray
        The array of values mapped to the domain [-1, 1].

    Examples
    --------
    >>> values = np.array([2, 2.25, 2.5, 2.75, 3])
    >>> a, b = 2, 3
    >>> inverse_linear_map(values, (a, b))
    array([-1. , -0.5,  0. ,  0.5,  1. ])
    """
    a, b = domain
    return 2 * (values - a) / (b - a) - 1


def domain_restricted_chebyt(n, domain):
    """
    Return the Chebyshev polynomial of the first kind of degree n, restricted to the domain [a, b].

    Parameters
    ----------
    n : integer
        The degree of the Chebyshev polynomial.
    domain : tuple
        The domain to restrict the Chebyshev polynomial to.

    Returns
    -------
    T : function
        The Chebyshev polynomial of the first kind of degree n, restricted to the domain [a, b].

    Examples
    --------
    >>> T = domain_restricted_chebyt(3, (-1, 1))
    >>> T(0)
    0.0
    """
    return lambda x: chebyt(n)(inverse_linear_map(x, domain))


def restircted_chebfit(f, degree=10, points=100, domain=(-1, 1)):
    """
    Return the Chebyshev coefficients of the best fit polynomial of degree degree to the function f, restricted to the domain domain.

    Parameters
    ----------
    f : function
        The input function.
    degree : integer, optional
        The degree of the best fit polynomial. The default is 10.
    points : integer, optional
        The number of points to use for the fit. The default is 100.
    domain : tuple, optional
        The domain to restrict the fit to. The default is (-1, 1).

    Returns
    -------
    coefficients : ndarray
        The Chebyshev coefficients of the best fit polynomial.

    Examples
    --------
    >>> None
    """
    x = np.linspace(-1, 1, points)
    x_mapped = linear_map(x, domain)
    y = f(x_mapped)

    return chebfit(x, y, degree - 1)
