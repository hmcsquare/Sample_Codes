from igraph import *
from math import ceil, floor
import random
from evaluation import *
import collections
import time

    
def sampling_FF_stat(graph, dist_original, indices, neighbor_dict, sample_size, pf, pb = 0,burn_in = 0, start_node = False):
    """
    Forest fire sampling with sampling ratio p, 
    forward burning prob pf, backward burning prob pb
    """ 
    if start_node == False:
        start_node = random.choice(indices)
        
    
    Q = collections.deque()
    V = collections.deque()
    L = set() #visited nodes
    count = 0
    while len(L) < sample_size:
        if count == 0:
            Q.append(start_node)
            V.append(start_node)
            v = Q.popleft()
        if Q:
            v = Q.popleft()
        if v not in L:
            L.add(v)
        neighbors = neighbor_dict[v]
        for u in neighbors:
            if u not in L:
                if random.random() < pf:
                        if u not in Q:
                            Q.append(u)
            if u in L:
                if random.random() < pb:
                    if u not in Q:
                        Q.append(u)

        count += 1
            
        
        
    V = list(L)
    V = V[burn_in:-1]
    samples = V
    degrees = [graph.degree(v) for v in V]
    weights = [1 / sample_size for v in V]
    avg_deg = mean(degrees)
    dist_est = deg_to_dist(degrees, weights)
    
    return TVD(dist_original, dist_est), NRMSE(dist_original, dist_est)

def sampling_RW_stat(graph, dist_original, indices, neighbor_dict, sample_size, burn_in = 0):
    """
    simple random walk sampling with re-weighting
    graph : igraph object to be sampled
    sample_size : sample size
    burn_in : burn in period with default = 0
    """
    
    samples = []
    weights = []
    start_node = random.choice(indices)
    current_node = start_node
    while len(samples) < sample_size + burn_in:
        neighbor = list(neighbor_dict[current_node])
        next_node = random.choice(neighbor)  # pick a random node from neighbor of current node
        current_node = next_node
        samples.append(current_node)
        weights.append(graph.degree(current_node) ** -1)
        
    samples = samples[burn_in:-1]  # discard burn-in samples
    weights = weights[burn_in:-1]
    degrees = [graph.degree(vertex) for vertex in samples] # degrees of samples
    Z = sum(weights)
    weights = [weight / Z for weight in weights] # normalizing weights
    avg_deg = 0
    for i in range(len(samples)):
        avg_deg += degrees[i] * weights[i]
    
    dist_est = deg_to_dist(degrees, weights)
    return TVD(dist_original, dist_est), NRMSE(dist_original, dist_est)

def sampling_MH_stat(graph, dist_original, indices, neighbor_dict, sample_size, burn_in = 0):
    """
    metropolis-hastings sampling
    graph : igraph object to be sampled
    sample_size : sample size
    burn_in : burn in period with default = 0
    """
    samples = []
    start_node = random.choice(indices)
    current_node = start_node
    
    while len(samples) < sample_size + burn_in:
        neighbor = list(neighbor_dict[current_node])
        next_node = random.choice(neighbor)
        kv = graph.degree(current_node) # degree of current node
        kw = graph.degree(next_node) # degree of next_node
        p = random.random()
        if p <= (kv / kw): # if accepted
            current_node = next_node
            
        else: # if rejected
            current_node = current_node
        samples.append(current_node)
            
    samples = samples[burn_in:-1] # discard burn-in samples
    degrees = [graph.degree(sample) for sample in samples]
    weights = [1/len(degrees) for i in samples]
    avg_deg = 0
    for i in range(len(samples)):
        avg_deg += degrees[i] * weights[i]
    dist_est = deg_to_dist(degrees, weights)
    return TVD(dist_original, dist_est), NRMSE(dist_original, dist_est)


def sampling_MHDA_stat(graph, dist_original, indices, neighbor_dict, sample_size, burn_in = 0):
    """
    metropolis hastings sampling with delayed acception
    graph : igraph object to be sampled
    sample_size : sample size
    burn_in : burn in period with default = 0
    """    
    samples = []
    start_node = random.choice(indices)
    current_node = start_node
    prev_node = start_node
    isFirst = True
    while len(samples) < sample_size + burn_in:
        neighbor = list(neighbor_dict[current_node])
        
        if isFirst: #in case of the first sample
            next_node = random.choice(neighbor)
            prev_node = current_node
            current_node = next_node
            isFirst = False
        else:
            next_node = random.choice(neighbor)
            ku = graph.degree(current_node) #degree of the current node 
            kw = graph.degree(prev_node) # degree of the previous node
            kv = graph.degree(next_node) # degree of the first candidate node
            p = random.random()
            if p < (ku / kv): # if the first candidate node is accepted
                if (next_node == prev_node) and (ku > 1): 
                    # if the node backtracked and the current degree is bigger than 1
                    temp = neighbor
                    temp.remove(prev_node) # N(current node) \ previous node 
                    k = random.choice(temp) # second candidate node
                    kk = graph.degree(k) # degree of the second candidate node
                    q = random.random()
                    if q <= min([1, min([1, (ku / kk) ** 2]) * max([1, (kv/ku) ** 2])]):
                        # if the second candidate is accepted
                        prev_node = current_node
                        current_node = k
                    else:
                        # if the second candidate is rejected --> in case of backtracking
                        prev_node = current_node
                        current_node = next_node
                else: # if the first candidate is not a backtracked node or degree of the current node is 1
                    prev_node = current_node
                    current_node = next_node
            else: # if the first candidate is rejected
                current_node = current_node
                prev_node = prev_node
        samples.append(current_node)
    samples = samples[burn_in:-1]
    degrees = [graph.degree(vertex) for vertex in samples]
    weights = [1/len(degrees) for i in samples]        
    avg_deg = 0
    for i in range(len(samples)):
        avg_deg += degrees[i] * weights[i]
    dist_est = deg_to_dist(degrees, weights)
    return TVD(dist_original, dist_est), NRMSE(dist_original, dist_est)
 

def sampling_HMC2_stat(graph, dist_original, indices, neighbor_dict, sample_size, block_size = 3, alpha_mu = 0.05, alpha_std = 0.02, burn_in_p = 0):
    
    samples = []
    weights = []
    
    s = floor(sample_size / block_size)
    
    start_node = random.choice(indices)
    current_node = start_node
    prev_node = start_node
    isFirst = True
    while len(samples) < sample_size:
        for i in range(block_size):
            alpha = 0
            while(alpha <= 0 or alpha > 1):
                alpha = random.normalvariate(alpha_mu, alpha_std)
            for j in range(s):
                burn_in = ceil(burn_in_p * s)
                neighbor = neighbor_dict[current_node]
                
                if isFirst:
                    next_node = random.choice(neighbor)
                    prev_node = current_node
                    current_node = next_node
                    isFirst = False
                else:
                    next_node = random.choice(neighbor)
                    ku = graph.degree(current_node)
                    kw = graph.degree(prev_node)
                    kv = graph.degree(next_node)
                    p = random.random()
                    if p < (ku / kv) ** alpha:
                        if (next_node == prev_node) and (ku > 1):
                            temp_neighbor = set(neighbor)
                            temp_neighbor.remove(prev_node)
                            k = random.choice(list(temp_neighbor))
                            kk = graph.degree(k)
                            q = random.random()
                            if q <= min(1, min(1, (ku / kk) ** 2) * max(1, (kv/ku) ** 2)) ** alpha:
                                prev_node = current_node
                                current_node = k
                            else:
                                prev_node = current_node
                                current_node = next_node
                            
                        else:
                            prev_node = current_node
                            current_node = next_node
                    else:
                        current_node = current_node
                        prev_node = prev_node
                if j > burn_in:
                    samples.append(current_node)
                    weights.append(graph.degree(current_node) ** (alpha - 1))
                if len(samples) >= sample_size:
                    break
            if len(samples) >= sample_size:
                break
    
    Z = sum(weights)
    weights = [weight / Z for weight in weights]
    degrees = [graph.degree(vertex) for vertex in samples]
    avg_deg = 0
    for i in range(len(samples)):
        avg_deg += degrees[i] * weights[i]
        
    dist_est = deg_to_dist(degrees, weights)
    return TVD(dist_original, dist_est), NRMSE(dist_original, dist_est)


if __name__ == "__main__":
    
    edgelist = sys.argv[1] + '_final'
    sample_size = int(sys.argv[2])
    g = get_graph(edgelist)
    print(g.summary())
    dist_original = deg_to_dist(g.degree())
    hashfile = sys.argv[1] + '_hash.p'
    indices, neighbor_dict = pickle.load(open(hashfile, 'rb'))
    TVD_FF, nrmse_FF = sampling_FF_stat(g, dist_original, indices, neighbor_dict, sample_size, 0.7) 
    TVD_RW, nrmse_RW = sampling_RW_stat(g, dist_original, indices, neighbor_dict, sample_size)
    TVD_MH, nrmse_MH = sampling_MH_stat(g, dist_original, indices, neighbor_dict, sample_size) 
    TVD_MHDA, nrmse_MHDA = sampling_MHDA_stat(g, dist_original, indices, neighbor_dict, sample_size)
    TVD_HMC2, nrmse_HMC2 = sampling_HMC2_stat(g, dist_original, indices, neighbor_dict, sample_size)
    
    print("TVD of forest fire sampling is : " + str(TVD_FF))
    print("TVD of random walk sampling is : " + str(TVD_RW))
    print("TVD of metropolis hastings sampling is : " + str(TVD_MH))
    print("TVD of metropolis hastings with delayed acceptance sampling is : " + str(TVD_MHDA))
    print("TVD of our method is : " + str(TVD_HMC2))
    print()
    print("nrmse of forest fire sampling is : " + str(nrmse_FF))
    print("nrmse of random walk sampling is : " + str(nrmse_RW))
    print("nrmse of metropolis hastings sampling is : " + str(nrmse_MH))
    print("nrmse of metropolis hastings with delayed acceptance sampling is : " + str(nrmse_MHDA))
    print("nrmse of our method is : " + str(nrmse_HMC2))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     