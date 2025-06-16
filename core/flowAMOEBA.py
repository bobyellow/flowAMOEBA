import time
import numpy as np
import random
from core.getFlowNeighbors import getFlowNeighborsContiguity
from core.spatstats import calculateGetisG

__all__ = ['execFlowAMOEBA']

def execFlowAMOEBA(AREAS1, AREAS2, FlowValue, NeiLvl,
                   significance=0.01, mc_reps=500, headtail_pct=0.1):
    """
    FlowAMOEBA: adaptation of AMOEBA for flow data (Tao & Thill 2019)

    Parameters:
    - AREAS1: origin areas (list of polygons)
    - AREAS2: destination areas (list of polygons)
    - FlowValue: dict of (O,D) -> [flow]
    - NeiLvl: spatial neighbor level code
    - significance: p-value threshold
    - mc_reps: number of Monte Carlo permutations
    - headtail_pct: fraction of flows to seed at each tail

    Returns:
    - multiline string: "O, D, ClusterID, FlowValue"
    """
    start = time.time()

    # 1) build flow neighborhood dictionary
    Wflow = getFlowNeighborsContiguity(AREAS1, AREAS2, FlowValue, NeiLvl)

    # 2) collect all keys and values
    keys = list(FlowValue.keys())
    vals = [FlowValue[k] for k in keys]
    N = len(keys)

    # 3) compute global stats on flows
    arr = np.array(vals, dtype=float)
    mean = np.nanmean(arr)
    std = np.nanstd(arr)

    # 4) head-tail seed selection
    sorted_keys = sorted(keys, key=lambda k: FlowValue[k])
    n_seed = max(1, int(headtail_pct * N))
    seeds = sorted_keys[:n_seed] + sorted_keys[-n_seed:]

    # 5) grow clusters
    generatedClusters = {}
    clusterG = {}
    for s in seeds:
        cluster = [s]
        Gcurr = calculateGetisG(cluster, mean, std, FlowValue, N)
        # iterative growth
        while True:
            prev = Gcurr
            # frontier neighbors
            frontier = set()
            for h in cluster:
                frontier.update(Wflow.get(h, []))
            frontier.difference_update(cluster)
            if not frontier:
                break
            # choose ordering
            if Gcurr >= 0:
                ordered = sorted(frontier, key=lambda x: FlowValue[x], reverse=True)
            else:
                ordered = sorted(frontier, key=lambda x: FlowValue[x])
            improved = False
            for nb in ordered:
                cand = cluster + [nb]
                Gnew = calculateGetisG(cand, mean, std, FlowValue, N)
                if (Gcurr >= 0 and Gnew > Gcurr) or (Gcurr < 0 and Gnew < Gcurr):
                    cluster, Gcurr = cand, Gnew
                    improved = True
                    break
            if not improved or Gcurr == prev:
                break
        generatedClusters[s] = cluster
        clusterG[s] = Gcurr

    # 6) rank seeds by |G| descent
    ranked = sorted(generatedClusters.keys(), key=lambda s: abs(clusterG[s]), reverse=True)

    # 7) Monte Carlo and assign cluster IDs
    output = {k: 0 for k in keys}
    pos_id, neg_id = 1, -1
    for s in ranked:
        cluster = generatedClusters[s]
        # skip if any member already assigned
        if any(output[h] != 0 for h in cluster):
            continue
        # permute and test
        better = 0
        for i in range(mc_reps):
            perm = random.sample(keys, N)
            idx = {keys[j]: perm[j] for j in range(N)}
            perm_cluster = [idx[h] for h in cluster]
            Gp = calculateGetisG(perm_cluster, mean, std, FlowValue, N)
            if (clusterG[s] >= 0 and Gp > clusterG[s]) or (clusterG[s] < 0 and Gp < clusterG[s]):
                better += 1
        p = better / float(mc_reps)
        if p <= significance:
            cid = pos_id if clusterG[s] > 0 else neg_id
            if clusterG[s] > 0:
                pos_id += 1
            else:
                neg_id -= 1
            # **suppress singleton clusters**
            if len(cluster) > 1:
                for h in cluster:
                    if output[h] == 0:
                        output[h] = cid

    elapsed = time.time() - start
    kept = list(output.values())
    pos_ids = sorted({cid for cid in kept if cid > 0})
    neg_ids = sorted({cid for cid in kept if cid < 0})
    print(f"flowAMOEBA completed in {elapsed:.2f}s")
    print(f"Retained {len(pos_ids)} hot clusters")
    print(f"Retained {len(neg_ids)} cold clusters")

    # 8) build output string
    lines = ["O, D, Cluster, Flow"]
    for k in keys:
        oid, did = k
        lines.append(f"{oid}, {did}, {output[k]}, {FlowValue[k]}")
    return "\n".join(lines)
