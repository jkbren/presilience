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
import warnings
import community


def modified_shannon_entropy(G, f, removal='random',
                             niter=50, return_stdv=False):
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
    for _ in range(niter):
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

    if return_stdv and niter > 4:

        return np.array(H_msh_mean).mean(), np.array(H_msh_mean).std()

    else:
        return np.array(H_msh_mean).mean()


def resilience(G, ntimes=2, rate=51, output_list=True, removal='random',
               H_std=True, niter=50):
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

    out_mean = []
    out_stdv = []
    for _ in range(ntimes):
        H_out_mean = []
        H_out_stdv = []
        for f in np.linspace(0, 1, rate):
            H_msh_mean, H_msh_stdv = modified_shannon_entropy(G, f, removal,
                                                              niter, H_std)
            H_out_mean.append(H_msh_mean)
            H_out_stdv.append(H_msh_stdv)
            # H_out.append(modified_shannon_entropy(G, f, removal))

        out_mean.append(np.array(H_out_mean))
        out_stdv.append(np.array(H_out_stdv))

    out_mean = np.array(out_mean).mean(axis=0)
    out_stdv = np.array(out_stdv).mean(axis=0)

    if output_list:
        return out_mean

    else:
        return 1 - sum(out_mean) / rate


def add_node(G, m, n, method='random', alpha=1.0):
    """
    Add node to the network according to the method supplied. This new node may
    be added randomly, preferentially based on degree, or using insights about
    empirical distributions of protein-specific interaction patterns (here
    named 'bio_smart', 'random', and 'degree').

    Params
    ------
    G (nx.Graph):
        the (protein-protein interaction) network in question.

    m (int):
        the number of edges that each new node brings to the network.

    method (str):
        the method of node-addition in question. can be any of the following:
            - 'random': adds node's edges randomly.
            - 'degree': adds edges preferentially based on degree.
            - 'bio_smart': adds edges based on gene_expression data.

    alpha (float):
        in case the method=='degree', alpha tunes the preferential attachment
        exponent, which makes node_i more / less likely to attaching its m
        edges to high-degree nodes.

    Returns
    -------
    G (nx.Graph):
        the graph with nodes added.

    """

    nodes = list(G.nodes())
    N = len(nodes)

    if method == 'random':
        # add random node
        probs = [1 / N for i in range(N)]
        eijs = np.random.choice(nodes, size=(m,), replace=False, p=probs)

        for node_j in eijs:
            G.add_edge(n, node_j)

        return G

    if method == 'degree':
        # add preferential attachment based on degree
        degrees = np.array(list(dict(G.degree()).values()))
        probs = (degrees**alpha) / sum(degrees**alpha)
        eijs = np.random.choice(nodes, size=(m,), replace=False, p=probs)

        for node_j in eijs:
            G.add_edge(n, node_j)

        return G

    if method == 'bio_smart':
        # add node preferentially by the node's gene expression attribute
        gex_dict = nx.get_node_attributes(G, 'gene_expression')
        gene_expression = np.array(list(gex_dict.values()))
        gene_exp_decile = np.percentile(gene_expression, 10)

        probs = (gene_expression) / sum(gene_expression)
        eijs = np.random.choice(nodes, size=(m,), replace=False, p=probs)

        for node_j in eijs:
            G.add_edge('added_protein_'+str(n), node_j)

        gene_exp_dict = nx.get_node_attributes(G, 'gene_expression')
        gene_exp_dict['added_protein_'+str(n)] = gene_exp_decile
        nx.set_node_attributes(G, gene_exp_dict, 'gene_expression')

        return G


def presilience(G, t=4, m=2, method='random', rate=100,
                ntimes=10, output_list=True, printt=True):
    """
    The 'presilience' is defined as the change in resilience (as calculated
    in Zitnik et al. (2019) using a modified Shannon entropy of the
    cluster size distribution in a network following uniform node removal)
    following the addition of a new node into the network. This new node
    may be added randomly, preferentially based on degree, or using
    insights about empirical distributions of protein-specific interaction
    patterns (named 'bio_smart', 'random', and 'degree').

    Parameters
    ----------
    G (nx.Graph):
        the protein-protein interaction network in question.

    t (int):
        the number of new nodes added (aka the number of timesteps) in the
        future that the presilience will be calculated.

    m (int):
        the number of edges that each new node brings to the network.

    method (str):
        the method of node-addition in question (can be either 'random'--adds
        node's edges randomly, 'degree'--adds edges preferentially based on
        degree, and 'bio_smart'-- which adds edges based on biological data).

    rate (int):
        the number of intervals between 0 and 1, which correspond to fractions
        of the network that are removed at each step.

    n_times (int):
        the number of runs that the algorithm goes through in order to arrive
        at the final (averaged) entropy value.

    output_list (bool):
        if True, returns list of resilience values. else return single number.

    Returns
    -------
    G (nx.Graph):
        the new graph with nodes added

    presilience (list or float):
        a list of resilience values or final resilience

    """

    Gx = G.copy()
    presilience = [resilience(Gx, ntimes, rate, output_list=False)]

    for new_node in range(t):
        if printt:
            print("\t Presilience t =", new_node)

        Gx = add_node(Gx, m, new_node, method)
        presilience.append(resilience(Gx, ntimes, rate, output_list=False))

    if output_list:
        return Gx, presilience
    else:
        return Gx, presilience[-1]


def presilience_mean(G, t=4, m=2, method='random', rate=40,
                     ntimes=10, output_list=True, n_iter=20, printt=True):
    """
    Runs the presilience algorithm several (n_iter) times.

    Parameters
    ----------
    G (nx.Graph):
        the protein-protein interaction network in question.

    t (int):
        the number of new nodes added (aka the number of timesteps) in the
        future that the presilience will be calculated.

    m (int):
        the number of edges that each new node brings to the network.

    method (str):
        the method of node-addition in question (can be either 'random'--adds
        node's edges randomly, 'degree'--adds edges preferentially based on
        degree, and 'bio_smart'-- which adds edges based on biological data).

    rate (int):
        the number of intervals between 0 and 1, which correspond to fractions
        of the network that are removed at each step.

    n_times (int):
        the number of runs that the algorithm goes through in order to arrive
        at the final (averaged) entropy value.

    output_list (bool):
        if True, returns list of resilience values. else, return single number.

    n_iter (int):
        number of iterations that go into creating the mean presilience.

    Returns
    -------
    presilience_mean (np.array):
        vector of resilience values of length t, averaged over n_iter times.

    """

    if output_list:
        pres = []
        for i in range(n_iter):
            if printt:
                print('Presilience run: %02i' % (i))

            _, pres_i = presilience(G, t, m, method, rate, ntimes, output_list)
            pres.append(pres_i)

        pres = np.array(pres)
        presilience_mean = pres.mean(axis=0)

    else:
        pres = []
        for i in range(n_iter):
            if printt:
                print('Presilience run: %02i' % (i))

            _, pres_i = presilience(G, t, m, method, rate, ntimes, output_list)
            pres.append(pres_i)

        pres = np.array(pres)
        presilience_mean = sum(pres) / len(pres)

    return presilience_mean


def modularience(G, t=4, m=2, method='random', output_list=True, printt=True):
    """
    The 'modularience' is defined as the change in modularity of the
    community-detected partition following the following the
    addition of a new node into the network. This new node may be
    added randomly, preferentially based on degree, or using insights
    about empirical distributions of protein-specific interaction
    patterns (named 'bio_smart', 'random', and 'degree').

    Parameters
    ----------
    G (nx.Graph):
        the protein-protein interaction network in question.

    t (int):
        the number of new nodes added (aka the number of timesteps) in the
        future that the presilience will be calculated.

    m (int):
        the number of edges that each new node brings to the network.

    method (str):
        the method of node-addition in question (can be either 'random'--adds
        node's edges randomly, 'degree'--adds edges preferentially based on
        degree, and 'bio_smart'-- which adds edges based on biological data).

    output_list (bool):
        if True, returns list of resilience values. else return single number.

    Returns
    -------
    G (nx.Graph):
        the new graph with nodes added.

    modularience (list or float):
        a list of modularity values or final modularity.

    """

    Gx = G.copy()
    partition = community.best_partition(Gx)
    modularit = [community.modularity(partition, Gx)]

    for new_node in range(t):
        if printt:
            print("\t Modularience t =", new_node)

        Gx = add_node(Gx, m, str(new_node), method)
        partition = community.best_partition(Gx)
        modularit.append(community.modularity(partition, Gx))

    if output_list:
        return Gx, modularit
    else:
        return Gx, modularit[-1]


def modularience_mean(G, t=4, m=2, method='random',
                      output_list=True, n_iter=20, printt=True):
    """
    Runs the modularience algorithm several (n_iter) times.

    Parameters
    ----------
    G (nx.Graph):
        the protein-protein interaction network in question.

    t (int):
        the number of new nodes added (aka the number of timesteps) in the
        future that the presilience will be calculated.

    m (int):
        the number of edges that each new node brings to the network.

    method (str):
        the method of node-addition in question (can be either 'random'--adds
        node's edges randomly, 'degree'--adds edges preferentially based on
        degree, and 'bio_smart'-- which adds edges based on biological data).

    output_list (bool):
        if True, returns list of resilience values. else return single number.

    Returns
    -------
    comm_mean (np.array):
        vector of modularity values of length t, averaged over n_iter times.

    """

    if output_list:
        comm = []
        for i in range(n_iter):
            if printt:
                print('Modularience run: %02i' % (i))

            _, comm_i = modularience(G, t, m, method, output_list)
            comm.append(comm_i)

        comm = np.array(comm)
        comm_mean = comm.mean(axis=0)

    else:
        comm = []
        for i in range(n_iter):
            if printt:
                print('Modularience run: %02i' % (i))

            _, comm_i = modularience(G, t, m, method, output_list)
            comm.append(comm_i)

        comm = np.array(comm)
        comm_mean = sum(comm) / len(comm)

    return comm_mean


def gene_expression_shuffle(G, p_permute):
    """
    Shuffles p_permute of the gene expression values. Each node in G
    should have a node attribute that corresponds to the gene expression
    of that protein. If we want to ask what the effect of *mere* values
    associated with each node (or the distribution of those values),
    this function becomes especially useful.

    Parameters
    ----------
    G (nx.Graph):
        the protein-protein interaction network in question.

    p_permute (float):
        the fraction of nodes in the network whose gene expression values
        will be shuffled among one another.

    Returns
    -------
    H (nx.Graph):
        the modified graph with new (shuffled) gene expression values.

    """

    gene_expression = nx.get_node_attributes(G, 'gene_expression')

    # permute the gene expression values among n nodes
    n_permute = int(p_permute * G.number_of_nodes())
    permutate_nodes = np.random.choice(list(G.nodes()), replace=False,
                                       size=n_permute)

    gene_expression_values = [gene_expression.get(i) for i in permutate_nodes]

    gene_expression_out = gene_expression.copy()
    np.random.shuffle(gene_expression_values)

    for p_i, perm in enumerate(permutate_nodes):
        gene_expression_out[perm] = gene_expression_values[p_i]

    H = G.copy()
    nx.set_node_attributes(H, gene_expression_out, 'gene_expression')

    return H
