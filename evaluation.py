from sampling import *
from utility import *
import pickle
import numpy as np

def deg_to_dist(degrees, weights = False):
    
    if weights == False:
        weights = [1 / len(degrees) for i in range(len(degrees))]
        dist = dict()
        for i in range(len(degrees)):
            if degrees[i] not in dist:
                dist[degrees[i]] = weights[i]
            else:
                dist[degrees[i]] += weights[i]

    else:
        assert len(degrees) == len(weights)
        dist = dict()
        for i in range(len(degrees)):
            if degrees[i] not in dist:
                dist[degrees[i]] = weights[i]
            else:
                dist[degrees[i]] += weights[i]
    return dist


def TVD(dist_original, dist_est):
    
    tvd = 0
    for i in dist_original:
        if i not in dist_est:
            px = 0
        else:
            px = dist_est[i]
        py = dist_original[i]
        tvd_i = abs(px - py)
        tvd += tvd_i
    tvd = tvd / 2
    
    return tvd
    
def NRMSE_(estimated_x, real_x):
    """
    calculate normalized rooted mean squared error
    """
    nrmse= np.sqrt((estimated_x -real_x)**2/real_x)
    
    return nrmse
 
def NRMSE_per_degree(dist_original, dist_est):
    
        
    nrmse_per_degree = dict()
    
    for i in dist_original:
        if i not in dist_est:
            px = 0
        else:
            px = dist_est[i]
        py = dist_original[i]
        
        nrmse = NRMSE_(px, py)
        
        nrmse_per_degree[i] = nrmse
        
    return nrmse_per_degree

def NRMSE(dist_original, dist_est):
    return mean(NRMSE_per_degree(dist_original, dist_est).values())
