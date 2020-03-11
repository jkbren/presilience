"""
utilities.py
--------------

Prospective resilience based off of [1], building on the resilience measure
that was introduced in [2]. This is used to simulate the change in resilience
following the addiiton of new nodes to a network, in this case ribosomal
protein-protein interaction networks.

[1] Klein, B., Holmér, L., Smith, K., Johnson, M., Swain, A., Stolp, L.,
Teufel, A., and Kleppe, A. (2020). Capturing the evolutionary capacity to
innovate via novel interactions in protein-protein interaction networks.

[2] Zitnik, M., Sosič, R., Feldman, M.W., and Leskovec, J. (2019). Evolution of
resilience in protein interactomes across the tree of life. Proceedings of the
National Academy of Sciences. 116, 10, 4426–4433. doi: 10.1073/pnas.1818013116.

author: Brennan Klein and Ludvig Holmér
email: brennanjamesklein@gmail.com
"""

import numpy as np


def softmax(A, k=1.0):
    """
    Calculates the softmax of a distribution, A,
    modulated by a precision term, k.

    Parameters
    ----------
    A (np.ndarray or list):
        the vector that you would like to softmax.

    k (float):
        when k is very negative, this function returns a
        uniform distribution. when high and positive, it
        resembles a delta funciton around A.max.

    Returns
    -------
    A (np.ndarray)
        the softmaxed version of the input vector.

    """
    A = np.array(A) if not isinstance(A, np.ndarray) else A
    A = A * k
    maxA = A.max()
    A = A - maxA
    A = np.exp(A)
    A = A / np.sum(A)

    return A


def average_every_n(xvec, yvec, n=2):
    """
    Utility function that spits out a smoothed x and y vector.

    Parameters
    ----------
    xvec, yvec (np.array):
        vectors of x and y data
    n (int):
        average every n terms together

    Returns
    -------
    out_x, out_y (np.array):
        two smoothed vectors according to however many n were specified

    """

    out_x = []
    out_y = []

    min_xdiff = xvec[1] - xvec[0]

    for i in range(0, len(xvec), n):
        xnumerat = 0
        ynumerat = 0

        if i + n <= len(xvec):
            for j in range(n):
                xnumerat += xvec[i+j]
                ynumerat += yvec[i+j]

            out_x.append(xnumerat / n)
            out_y.append(ynumerat / n)

    out_x = np.array(out_x)
    out_y = np.array(out_y)

    if n > 1:
        out_x = out_x + min_xdiff / n

    return out_x, out_y


def get_nodesizes(G, attach_type, ns=30):
    """
    Simple function that returns a vector of node sizes, based on some
    attribute, for pretty plotting.

    Parameters
    ----------
    G (nx.Graph):
        the network to be plotted.

    attach_type (str):
        string corresponding to attribute that should be used for sizing nodes.

    ns (float):
        baseline size of nodes, initialized to ns=30.

    Returns
    -------
    (np.ndarray):
        array of node sizes, indexed the same way as G.nodes() is.

    """

    if attach_type == 'none':
        return np.array([ns]*G.number_of_nodes())

    if attach_type == 'random':
        return np.random.uniform(0.5*ns, 1.5*ns, G.number_of_nodes())

    if attach_type == 'gene-expression':
        return ns * np.random.lognormal(0, 0.85, G.number_of_nodes())

    degs = np.array(list(dict(G.degree()).values()))
    if attach_type == 'degree':
        return 2 + np.array(list(dict(G.degree()).values()))**1.75

    if attach_type == 'inverse-degree':
        return (-degs-min(-degs)+1)**1.5
