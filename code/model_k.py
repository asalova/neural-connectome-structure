import numpy as np
import pandas
import sys
import copy
from scipy.spatial.distance import squareform

# Input: organism and constraint
name = sys.argv[1] # 'fly', 'mouse', or 'human'
constraint_name = sys.argv[2] # 'all' or 'cont'

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

def degree_preserving_rand(A_c, deg_dist, is_conn, tol=1e-2, max_iter=1000, check_timestep=100):
    '''
    Performs degree preserving randomization with constraints.

    Args:
        A_c (np array): Constraint matrix.
        deg_dist (np array): Degree distribution of the connectome.
        is_conn (np array): Vector indicating whether the node is connected.
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
    deg_dist = deg_dist[is_conn]
    size = len(A_c)
        
    alphas_opt = lagrange_mult_degree_preserving_rand(deg_dist, A_c, size, tol, max_iter, check_timestep)
    alphas = np.zeros(n)
    alphas[is_conn] = alphas_opt[:,0]

    alphas = np.expand_dims(alphas, 1)
    probs_mat = A_cont/(1 + alphas @ alphas.T)
    probs_mat[not_conn] = 0
    probs_mat[:, not_conn] = 0
                
    return alphas, probs_arr       

def lagrange_mult_degree_preserving_rand(deg_dist, A, size, tol, max_iter, check_timestep):
    '''
    Calculates the Lagrange multipliers.
    
    Returns:
        np array: Lagrange multipliers of the maximum entropy model.
    '''

    alphas_old = np.ones((size,1))
    for i in range(max_iter):
        if i%check_timestep==0:
            print(f'{i} iterations', end="\r", flush=True)
        alphas_new = (1/deg_dist) * np.sum(A/(alphas_old.T+1/alphas_old), axis=1)
        if i%check_timestep==0:
            deg_dist_new = np.sum(A/(1+alphas_old*alphas_old.T), axis=1)
            max_err = np.max(np.abs(deg_dist - deg_dist_new))
            if max_err<tol:
                print(f'Completed after {i+1} iterations, maximum error {max_err}')
                break
        alphas_old = np.expand_dims(alphas_new, 1)
    if i==(max_iter-1):
        print(f'Stopped after {i+1} iterations, maximum error {max_err}')
    return alphas_old

degree_preserving_rand(A_constr, deg_dist, is_conn)