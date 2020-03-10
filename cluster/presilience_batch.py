import json
import random
import sys

import networkx as nx
import numpy as np
import datetime as dt
from presilience import modularience_mean,\
                        gene_expression_shuffle,\
                        presilience_mean

# Start simulation
args = list(sys.argv)
print(args)
# use index to import network
network_index = str(args[1])
if network_index == 'sce':
    G = nx.read_graphml('data/G_sce.graphml')

if network_index == 'eco':
    G = nx.read_graphml('data/G_eco.graphml')

if network_index == 'hsa':
    G = nx.read_graphml('data/G_hsa.graphml')

# simulation parameterizations
links_per_new = int(args[2])
method_used = str(args[3])
iterations = int(args[4])
seed_num = int(args[5])

# set seed
random.seed(seed_num)

timesteps_out = 20
noise_interval = np.linspace(0, 1, 11).round(3)
r = 40
nrep_presi = 4

# these will be the outputs
pres = {}
comm = {}

print("Starting at", dt.datetime.now())

if method_used == 'bio_smart':
    for i in noise_interval:
        if i > 0:
            H = gene_expression_shuffle(G, i)

        else:
            H = G.copy()

        pres_h = presilience_mean(H, t=timesteps_out, m=links_per_new,
                                  method=method_used, rate=r,
                                  ntimes=nrep_presi, output_list=True,
                                  n_iter=iterations)
        pres[i] = list(pres_h)
        # pres.append(list(pres_h))

        comm_h = modularience_mean(H, t=timesteps_out, m=links_per_new,
                                   method=method_used, output_list=True,
                                   n_iter=iterations)
        comm[i] = list(comm_h)
        # comm.append(list(comm_h))


else:
    H = G.copy()

    pres_h = presilience_mean(H, t=timesteps_out, m=links_per_new,
                              method=method_used, rate=r,
                              ntimes=nrep_presi, output_list=True,
                              n_iter=iterations)

    comm_h = modularience_mean(H, t=timesteps_out, m=links_per_new,
                               method=method_used, output_list=True,
                               n_iter=iterations)

    pres[noise_interval[0]] = list(pres_h)
    # pres.append(list(pres_h))
    comm[noise_interval[0]] = list(comm_h)
    # comm.append(list(comm_h))

data_m = {links_per_new: {method_used: {'pres': pres, 'comm': comm}}}
out_dict = {network_index: data_m}

print("Done at", dt.datetime.now())

json = json.dumps(out_dict, indent=4)
f = open("out/presil_%s_%s_%02i_%02i_%05i.json" % (network_index,
                                                   method_used,
                                                   links_per_new,
                                                   iterations,
                                                   seed_num), "w")

f.write(json)
f.close()
