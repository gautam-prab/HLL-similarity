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
        hll.insert(seq)
    return hll

folder_path = 'Data/0.5x'
sketches = []
species = []
i = 0
for filename in glob.glob(os.path.join(folder_path, '*.fasta')):
    with open(filename) as myFile:
        inputLines = myFile.readlines()
    species.append(filename)
    totalReads = int(len(inputLines) / 2)
    sketches.append( get_HLL(inputLines, totalReads))
    i += 1

similarity_matrix = np.zeros((10,10))

print(sketches)

for i in range(10):
    for j in range(i, 10):
        h1 = sketches[i]
        h2 = sketches[j]
        union = Similarity.union(h1, h2).cardinality()
        intersection = Similarity.intersection(h1, h2)
        jaccard = intersection / union
        similarity_matrix[i, j] = jaccard

np.set_printoptions(precision = 2, linewidth = 200)
print('Similarity Matrix:')
print(similarity_matrix)

rankings = GenomeRankings.rank_genomes(similarity_matrix)
print('Rankings:')
for i in range(rankings.shape[0]):
    print('{}\t\t{}'.format(species[i],str(rankings[i,:])))

# ground_truth = np.array()
# acc = GenomeRankings.compare_rankings(rankings, ground_truth)
