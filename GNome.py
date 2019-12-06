import Cardinality
import os
import numpy as np
from HLL import HLL
import Similarity
import glob
import GenomeRankings


def get_HLL(reads, numReads):
    hll = HLL(12)
    for i in range(numReads):
        seq = reads[(i*2)+1].rstrip().upper()
        # add all 25-mers to HLL
        for i in range(0,len(seq)-25):
            hll.insert(seq[i:i+25])
    return hll

folder_path = 'Data/50x'
sketches = []
species = []
i = 0

print('Reading Files...')
files = glob.glob(os.path.join(folder_path, '*.fasta'))
files.sort()
for filename in files:
    with open(filename) as myFile:
        inputLines = myFile.readlines()
    species.append(filename)
    print('\tReading file: {}'.format(filename))
    totalReads = int(len(inputLines) / 2)
    sketches.append( get_HLL(inputLines, totalReads))
    i += 1

jaccard_matrix = np.zeros((10,10))
sd_matrix = np.zeros((10,10))
forbes_matrix = np.zeros((10,10))

print('Calculating pairwise similarities...')
for i in range(10):
    for j in range(i, 10):
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

# np.set_printoptions(precision = 2, linewidth = 200)
# print('Similarity Matrix:')
# print(similarity_matrix)

jaccard_rankings = GenomeRankings.rank_genomes(jaccard_matrix)
sd_rankings = GenomeRankings.rank_genomes(sd_matrix)
forbes_rankings = GenomeRankings.rank_genomes(forbes_matrix)

# print('Rankings:')
# for i in range(rankings.shape[0]):
#     print('{}\t\t{}'.format(species[i],str(rankings[i,:])))

ground_truth = np.array([[0., 9., 4., 2., 1., 9., 9., 9., 9., 3.],
       [9., 0., 9., 9., 9., 9., 9., 9., 9., 9.],
       [3., 9., 0., 2., 1., 9., 9., 9., 9., 9.],
       [2., 9., 3., 0., 1., 9., 9., 9., 9., 4.],
       [2., 9., 3., 1., 0., 9., 9., 9., 9., 4.],
       [9., 9., 9., 9., 9., 0., 1., 9., 9., 9.],
       [9., 9., 9., 9., 9., 1., 0., 9., 9., 9.],
       [9., 9., 9., 9., 9., 9., 9., 0., 1., 9.],
       [9., 9., 9., 9., 9., 9., 9., 1., 0., 9.],
       [2., 9., 9., 3., 1., 9., 9., 9., 9., 0.]]) # from ANI data

jaccard_acc = GenomeRankings.compare_rankings(jaccard_rankings, ground_truth)
sd_acc = GenomeRankings.compare_rankings(sd_rankings, ground_truth)
forbes_acc = GenomeRankings.compare_rankings(forbes_rankings, ground_truth)
print('Jaccard Similarity Accuracy: {}'.format(jaccard_acc))
print('Sorenson-Dice Similarity Accuracy: {}'.format(sd_acc))
print('Forbes Similarity Accuracy: {}'.format(forbes_acc))

import pickle
with open('50x.sketch', 'wb') as sketchfile:
  pickle.dump(sketches, sketchfile)
