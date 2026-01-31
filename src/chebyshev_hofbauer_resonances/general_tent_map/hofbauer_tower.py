import numpy as np


def intersect(domain1, domain2):
    """
    Compute the intersection of two domains.

    Parameters
    ----------
    domain1 : tuple
        Tuple of (start, end) for the first domain.
    domain2 : tuple
        Tuple of (start, end) for the second domain.

    Returns
    -------
    tuple or None
        Tuple of (start, end) for the intersected domain, or None if the domains do not overlap.
    """
    start = max(domain1[0], domain2[0])
    end = min(domain1[1], domain2[1])
    if start < end:
        return (start, end)
    else:
        return None


def construct_adj_matrix(domains, ranges):
    """
    Construct the adjacency matrix from the domains and ranges.

    Parameters
    ----------
    domains : list
        List of tuples of the domains for each function.
    ranges : list
        List of ranges. Each range is either a tuple or None.

    Returns
    -------
    ndarray
        Adjacency matrix where entry (i, j) is 1 if the range of function j equals domain i.
    """
    n = len(domains)
    adj_matrix = np.zeros((n, n), dtype=int)

    for j, domain in enumerate(domains):
        range = ranges[j]
        if range is not None:
            for i, target_domain in enumerate(domains):
                if target_domain == range:
                    adj_matrix[i, j] = 1

    return adj_matrix


def create_adjacency_matricies(function_domains, functions, depth=1):
    """
    Compute the Hofbauer tower adjacency matrices for a piecewise map.

    Iteratively constructs the Hofbauer tower by tracking how domains map through
    each branch of the piecewise function. Starting from the unit interval [0, 1],
    this function intersects working domains with each function's domain, computes
    the resulting ranges, and builds adjacency matrices that encode the transitions
    between domains for each function branch.

    Parameters
    ----------
    function_domains : list
        List of tuples (start, end) specifying the domain for each function branch.
    functions : list
        List of callable functions, one for each domain.
    depth : int, optional
        Number of iterations to build the tower. Default is 1.

    Returns
    -------
    complete_domains : list
        List of all domain tuples discovered during tower construction.
    adj_matrices : list
        List of adjacency matrices, one for each function branch. Each matrix
        has shape (n, n) where n is the number of complete domains.
    """

    complete_domains = []
    working_domains = [(0, 1)]
    new_domains = []
    function_ranges = {i: [] for i in range(len(functions))}

    for _ in range(depth):
        while working_domains:
            current_domain = working_domains.pop()
            if current_domain in complete_domains:
                continue
            for i, (func_domain, func) in enumerate(zip(function_domains, functions)):
                intersected_domain = intersect(current_domain, func_domain)
                if intersected_domain is not None:
                    range_start = func(intersected_domain[0])
                    range_end = func(intersected_domain[1])
                    new_range = (
                        min(range_start, range_end),
                        max(range_start, range_end),
                    )
                    new_domains.append(new_range)
                    function_ranges[i].append(new_range)
                else:
                    function_ranges[i].append(None)

            complete_domains.append(current_domain)
        working_domains = new_domains
        new_domains = []

    adj_matrices = [
        construct_adj_matrix(complete_domains, function_ranges[i])
        for i in range(len(functions))
    ]

    return complete_domains, adj_matrices
