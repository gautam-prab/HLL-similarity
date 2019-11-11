from HLL import HLL
import Similarity
import Random_gen

import numpy as np
import matplotlib.pyplot as plt

def main():
    h1 = HLL(12)
    h2 = HLL(12)

    num_jacs = 100
    jaccards = np.linspace(0,0.5,num_jacs)
    plot = np.zeros(num_jacs)
    for i in range(num_jacs):
        jac = jaccards[i]
        # do 1 trials for each Jaccard
        trials = 100
        results = np.zeros(trials)
        for j in range(trials):
            # generate 10,000 reads for each
            a,b,exp_jaccard = Random_gen.generate_reads(jac, 10000, 10000, 50)

            for r in a:
                h1.insert(r)
            for r in b:
                h2.insert(r)

            union = Similarity.union(h1,h2).cardinality()
            intersection = Similarity.intersection_inclusion_exclusion(h1,h2)

            obs_jaccard = intersection/union
            error = (obs_jaccard - exp_jaccard)
            results[j] = error

        plot[i] = np.mean(results)

    plt.scatter(jaccards,plot)
    plt.show()

if (__name__ == '__main__'):
    main()
