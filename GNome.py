"""
GNome.py: Similarity estimation for simulated bacterial genome reads

Usage: accepts one argument which is the path to read Data
       python GNome.py Data/0.5x
"""

import Cardinality
import os
import sys
import numpy as np
from HLL import HLL
import Similarity
import glob
import GenomeRankings
import pickle

"""
get_HLL():
Generates and returns HLL by inserting 25-mers from reads

Input: reads from FASTA and number of reads in file
Returns: filled HLL from set of reads
"""
def get_HLL(reads, numReads):
    hll = HLL(12)
    for i in range(numReads):
        seq = reads[(i*2)+1].rstrip().upper()
        # add all 25-mers to HLL
        for i in range(0,len(seq)-25):
            hll.insert(seq[i:i+25])
    return hll

def get_sketches(folder_path):
    sketches = []
    species = []
    i = 0

    files = glob.glob(os.path.join(folder_path, '*.fasta'))
    files.sort()
    for filename in files:
        with open(filename) as myFile:
            inputLines = myFile.readlines()
        species.append(filename)
        print('\tReading file: {}'.format(filename))
        totalReads = int(len(inputLines) / 2)
        sketches.append(get_HLL(inputLines, totalReads))
        i += 1

    return sketches, species

"""
calculate_rankings():
Creates similarity matrices for three separate similarity coefficients
and calls rank_genomes to create vectors based on similarities

Input: array of sketches
Returns: matrix of each genomes similarity rankings
"""
def calculate_rankings(sketches):
    n_genomes = len(sketches)

    jaccard_matrix = np.zeros((n_genomes,n_genomes))
    sd_matrix = np.zeros((n_genomes,n_genomes))
    forbes_matrix = np.zeros((n_genomes,n_genomes))

    print('Calculating pairwise similarities...')
    for i in range(n_genomes):
        for j in range(i, n_genomes):
            h1 = sketches[i]
            h2 = sketches[j]
            union = Similarity.union(h1, h2).cardinality()
            a_excl, b_excl, intersection = Similarity.getJointEstimators(h1, h2)
            a = h1.cardinality()
            b = h2.cardinality()

            jaccard = intersection / union
            sd = 2*intersection/(a + b)
            forbes = (intersection*union)/(intersection*union + 1.5*a_excl*b_excl)

            jaccard_matrix[i, j] = jaccard
            sd_matrix[i, j] = sd
            forbes_matrix[i, j] = forbes

    jaccard_rankings = GenomeRankings.rank_genomes(jaccard_matrix)
    sd_rankings = GenomeRankings.rank_genomes(sd_matrix)
    forbes_rankings = GenomeRankings.rank_genomes(forbes_matrix)
    return jaccard_rankings, forbes_rankings, sd_rankings

"""
get_ground_truth():

Returns matrix of ground_truth similarities generated from enviomics ANI
"""
def get_ground_truth():
    return np.array([[0., 9., 4., 9., 2., 1., 9., 9., 9., 3.],
       [9., 0., 9., 9., 9., 9., 9., 9., 9., 9.],
       [3., 9., 0., 9., 2., 1., 9., 9., 9., 9.],
       [9., 9., 9., 0., 9., 9., 1., 9., 9., 9.],
       [2., 9., 3., 9., 0., 1., 9., 9., 9., 4.],
       [2., 9., 3., 9., 1., 0., 9., 9., 9., 4.],
       [9., 9., 9., 1., 9., 9., 0., 9., 9., 9.],
       [9., 9., 9., 9., 9., 9., 9., 0., 1., 9.],
       [9., 9., 9., 9., 9., 9., 9., 1., 0., 9.],
       [2., 9., 9., 9., 3., 1., 9., 9., 9., 0.]]) # from ANI data

"""
For each genome in given file, generates k-mers from reads and sketches HLL.
Then, calculates similarity ranking matrices for Jaccard, Forbes, and
Sorenson-Dice and compares to ground truth matrix.

Output: Similarity accuracies (percent matches between similarity matrix and
ground truth)
"""
def main():
    if sys.argv[1].rstrip()[-7:] == '.sketch':
        sketches = pickle.load(open(sys.argv[1], 'rb'))

    else:
        folder_path = sys.argv[1]  # path: 'Data/0.5x' or 5x/50x for other coverages
        print('Reading Files...')
        sketches, species = get_sketches(folder_path)


    jaccard_rankings, forbes_rankings, sd_rankings = calculate_rankings(sketches)
    ground_truth = get_ground_truth()

    jaccard_acc = GenomeRankings.compare_rankings(jaccard_rankings, ground_truth)
    sd_acc = GenomeRankings.compare_rankings(sd_rankings, ground_truth)
    forbes_acc = GenomeRankings.compare_rankings(forbes_rankings, ground_truth)
    print('Jaccard Similarity Accuracy: {}'.format(jaccard_acc))
    print('Sorenson-Dice Similarity Accuracy: {}'.format(sd_acc))
    print('Forbes Similarity Accuracy: {}'.format(forbes_acc))

    outfile = input('Output Sketch Filename (enter for none): ')
    if outfile.rstrip() != '':
        with open(outfile+'.sketch', 'wb') as sketchfile:
          pickle.dump(sketches, sketchfile)

if (__name__ == '__main__'):
    main()
