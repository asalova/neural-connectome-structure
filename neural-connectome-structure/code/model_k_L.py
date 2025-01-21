import numpy as np
import pandas
import sys
import copy
import os
from os.path import exists
from scipy.spatial.distance import squareform, pdist

# Input: organism and constraint
name = sys.argv[1] # 'fly', 'mouse', or 'human'
constraint_name = sys.argv[2] # 'all' or 'cont'

# Soma sizes are recorded for each organism
soma_size = {}
soma_size['fly'] = 2470
soma_size['mouse'] = 5340
soma_size['human_all'] = 7600

# Distance exponent parameters that match the total distance
exp_params = {}
exp_params['fly'] = 9.055
exp_params['mouse'] = 8.716
exp_params['human_all'] = 9.985
exp_param = exp_params[name]

# Loading neuron-level information
somas = pandas.read_cs­v(f'{name}_data/{name}_extended_neuron_info.csv')
# Loading the undirected unweighted synaptic network (connectome)
syn_undir = pandas.read_csv(f'{name}_data/{name}_weighted_connectome.csv')

# Getting degree distribution
deg_dist = somas.degree
is_conn = deg_dist>0

# Converting edge list to matrix 
def undir_mat_from_df(net_df, n):
    A = np.zeros((n,n))
    A[net_df.i, net_df.j] = 1
    A = A + A.T
    return A

# Loading the constraints. Even if there's no edge constraint ('all'), 
# self links need to be excluded.
if constraint_name=='cont':
    constr = pandas.read_csv(f'{name}/{name}_contactome.csv')
    A_constr = undir_mat_from_df(constr, len(somas))
if constraint_name == 'all':
    A_constr = - np.eye(len(somas)) + 1

# Calculating pairwise distances in soma size units.
pos_vect = np.array(somas[['x_cm', 'y_cm', 'z_cm']])
D = pdist(pos_vect)/soma_size[name]

def degree_preserving_rand(A_c, deg_dist, is_conn, D, exp_param, tol=1e-2, max_iter=1000, check_timestep=100):
    '''
    Performs degree preserving randomization with distance decay
    and hard constraints.

    Args:
        A_c (np array): Constraint matrix.
        deg_dist (np array): Degree distribution of the connectome.
        is_conn (np array): Vector indicating whether the node is connected.
        D (np array): Pairwise distances between neurons.
        exp_param (float): exponential decay parameter in max ent model.
        tol (float): Maximum allowed error in degree of the model.
        max_iter (int): Maximum number of iterations.
        check_timestep (int): Number of steps to print updates.
    
    Returns:
        alphas (np array): Lagrange multipliers of the max ent model.
        probs_arr (np array): Flattened edge probabilities.
    '''

    A_cont = copy.copy(A_c)
    not_conn = np.invert(is_conn)
    A_cont[not_conn] = 0
    A_cont[:, not_conn] = 0
 
    # Preserve degree distribution
    n = len(A_c)
    A_c = A_c[is_conn][:, is_conn]
    D_c = squareform(D)[is_conn][:, is_conn]
    deg_dist = deg_dist[is_conn]
    size = len(A_c)
    
    D_exp = np.exp(exp_param * D_c)
        
    alphas_opt = lagrange_mult_degree_preserving_rand(deg_dist, A_c, size, D_exp, tol, max_iter, check_timestep)
    alphas = np.zeros(n)
    alphas[is_conn] = alphas_opt[:,0]

    alphas = np.expand_dims(alphas, 1)
    probs_mat = A_cont/(1 + np.exp(exp_param * squareform(D)) * (alphas @ alphas.T))
    probs_mat[not_conn] = 0
    probs_mat[:, not_conn] = 0
    probs_arr = squareform(probs_mat)
                
    return alphas, probs_arr

def lagrange_mult_degree_preserving_rand(deg_dist, A, size, D_exp, tol, max_iter, check_timestep=10):
    '''
    Calculates the Lagrange multipliers.
    
    Returns:
        np array: Lagrange multipliers of the maximum entropy model.
    '''

    alphas_old = np.ones((size,1))
    for i in range(max_iter):
        if i%check_timestep==0:
            print(f'{i} iterations', end="\r", flush=True)
        #print(np.shape(alphas_old.T * D_exp))
        alphas_new = (1/deg_dist) * np.sum(A/(alphas_old.T * D_exp + 1/alphas_old), axis=1)
        if i%check_timestep==0:
            deg_dist_new = np.sum(A/(1+alphas_old * alphas_old.T * D_exp), axis=1)
            max_err = np.max(np.abs(deg_dist - deg_dist_new))
            if max_err<tol:
                print(f'Completed after {i+1} iterations, maximum error {max_err}', flush=True)
                break
            else:
                print(f'Err = {max_err}', flush=True)
        alphas_old = np.expand_dims(alphas_new, 1)
        if np.isnan(max_err):
            break
    if i==(max_iter-1):
        print(f'Stopped after {i+1} iterations, maximum error {max_err}', flush=True)
    return alphas_old

degree_preserving_rand(A_constr, deg_dist, is_conn, D, 1/exp_param)
