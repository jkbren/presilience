"""
presilience.py
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

import networkx as nx
import numpy as np
from scipy import stats
import scipy as sp
import warnings


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


def resilience(G, ntimes=10, rate=50, output_list=True, removal='random'):
    """
    The resilience of a network, G, is defined as the Shannon
    entropy of the cluster size distribution of a graph at a given
    time, t. By repeatedly removing a fraction of random nodes at each
    timestep, we observe the change of this distribution and, as such
    the entropy, H_sh. The resilience is calculated as the 1-sum(H_sh).

    Parameters
    ----------
    G (nx.Graph):
        the graph in question.

    n_times (int):
        the number of runs that the algorithm goes through in order to arrive
        at the final (averaged) entropy value.

    rate (int):
        the number of intervals between 0 and 1, which correspond to fractions
        of the network that are removed at each step.

    output_list (bool):
        if True, returns a list of resilience values. else, returns one value.
.
    removal (str):
        method of node-removal. for now this only includes 'random', but one
        can imagine a number of ways to systematically bias the node-removal
        process (e.g. based preferentially on degree, etc.)

    Returns
    -------
    out_mean (list or float):
        if output_list==True, this function returns a list of length = rate
        entropy values, which forms the curve that is used to calculate the
        network resilience. else, it returns 1 - sum(out_mean)/rate.

    """

    out = []
    for _ in range(ntimes):
        H_out = []
        for f in np.linspace(0, 1, rate):
            H_out.append(modified_shannon_entropy(G, f, removal))

        out.append(np.array(H_out))

    out = np.array(out)
    out_mean = out.mean(axis=0)

    if output_list:
        return out_mean

    else:
        return 1 - sum(out_mean)/rate


def modified_shannon_entropy(G, f, removal='random',
                             ntimes=10, return_stdv=False):
    """
    After a fraction, f, nodes have been 'removed' from the network (i.e., they
    have become disconnected isolates, such that the number of nodes does not
    change), this measure corresponds to the entropy of the distribution of
    component sizes, P(C). That is, P(c_i) corresponds to the probability that
    a randomly selected node will be in component c_i; there is a higher
    probability that a randomly-selected node will be in larger components.

    Parameters
    ----------
    G (nx.Graph):
        the graph in question.

    f (float):
        the fraction of N nodes in the network to be removed.

    removal (str):
        the method by which nodes are removed from the network.

    ntimes (int):
        the number of times to remove nodes to compute the mean entropy.

    return_stdv (bool):
        if True, this function returns a mean entropy as well as the standard
        deviation of the entropy after running ntimes. if False, this function
        just returns the mean value for the entropy.

    Returns
    -------
    H_msh_mean (float):
        the modified Shannon entropy averaged over ntimes runs. if return_stdv,
        this function also returns H_msh_stdv.

    """
    H_msh_mean = []
    for _ in range(ntimes):
        out_H = []
        N = G.number_of_nodes()
        leading_term = -1 / np.log2(N)
        p_i_unif = 1 / N

        G_f = G.copy()

        if removal == 'random':
            remove_nodes = [i for i in G_f.nodes() if np.random.rand() < f]

        else:
            warnings.warn("Only implemented for *random*. switching to that.")
            remove_nodes = [i for i in G_f.nodes() if np.random.rand() < f]

        curr = 0

        if len(remove_nodes):
            for node_r in remove_nodes:
                G_f.remove_node(node_r)
                curr += p_i_unif * np.log2(p_i_unif)

        Cs = list(nx.connected_components(G_f))
        for Ci in Cs:
            p_i = len(Ci) / N
            curr += p_i * np.log2(p_i)

        out_H.append(curr)

        H_msh = np.abs(leading_term * np.sum(out_H))
        H_msh_mean.append(H_msh)

    H_msh_mean = np.array(H_msh_mean).mean()

    if return_stdv and ntimes > 4:
        H_msh_stdv = np.array(H_msh_mean).mean()
        return H_msh_mean, H_msh_stdv

    else:
        return H_msh_mean
