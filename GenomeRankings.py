"""
GenomeRanking.py: Ranking genomes and calculating correlations from ranking matrices

This file supplies helper methods for GNome.py
"""
import numpy as np
from scipy.stats import rankdata
from scipy.stats import kendalltau

"""
rank_genomes():
Rank genomes based on similarity
Input is a square upper triangular matrix of similarity values of a similar form:
    genome 1    1   s(1,2)  s(1,3)
    genome 2    0   1       s(2,3)
    genome 3    0   0       1
(such as one from GNome.calculate_rankings())

Output is a matrix of each genome's rankings:
    genome 1    rank(1) rank(2) rank(3)
    genome 2    rank(1) rank(2) rank(3)
    genome 3    rank(1) rank(2) rank(3)
"""
def rank_genomes(sim_matrix):
    # copy upper triangular matrix into both triangles for ease of use
    sim_matrix = sim_matrix + sim_matrix.T - np.diagflat(sim_matrix.diagonal())

    # pairwise-rank each genome based on similarities
    # lower ranks are better matches
    # each row is a ranking vector for that corresponding genome
    # ties are ranked with the minimum ranking

    rankings = np.zeros(sim_matrix.shape)
    # iterate through rows since scipy doesn't have basic matrix functionality
    for idx,row in enumerate(sim_matrix):
        rankings[idx,:] = len(row) - rankdata(row, method='min')
    return rankings

"""
compare_rankings():

Compare two ranking matrices using average Kendall's tau distance

Input: two ranking matrices from the rank_genomes() function
Output: a float between -1 and 1; 1 indicates perfect agreement, -1 indicates perfect disagreement"""
def compare_rankings(mat1, mat2):
    n = mat1.shape[0]
    kt_sum = 0
    for idx in range(n):
        tau,_ = kendalltau(mat1[idx,:], mat2[idx,:])
        kt_sum += tau
    return kt_sum / n
