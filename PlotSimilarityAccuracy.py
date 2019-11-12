from HLL import HLL
import Similarity
import Random_gen

import sys

import numpy as np
import matplotlib.pyplot as plt

"""USAGE:
    python PlotSimilarityAccuracy.py jac to modulate jaccard
    python PlotSimilarityAccuracy.py card to modulate cardinality
"""

def main(arg):

    if arg == 'jac':
        num_jacs = 50
        jaccards = np.linspace(0,0.2,num_jacs)
        plot = np.zeros(num_jacs)
        for i in range(num_jacs):
            print('Starting '+str(i+1)+' out of '+str(num_jacs))
            jac = jaccards[i]
            # do 1 trials for each Jaccard
            trials = 5
            results = np.zeros(trials)
            for j in range(trials):
                # generate 10,000 reads for each
                a,b,exp_jaccard = Random_gen.generate_reads(jac, 10000, 10000, 40)
                h1 = HLL(12)
                h2 = HLL(12)

                for r in a:
                    h1.insert(r)
                for r in b:
                    h2.insert(r)

                union = Similarity.union(h1,h2).cardinality()
                intersection = Similarity.intersection_inclusion_exclusion(h1,h2)

                obs_jaccard = intersection/union
                error = 100*(obs_jaccard - exp_jaccard)/(1+exp_jaccard)
                results[j] = error

            plot[i] = np.mean(results)

        print(plot)
        plt.title('Jaccard Accuracy for 10,000 reads of length 40')
        plt.xlabel('Jaccard Index')
        plt.ylabel('% Error (mean of 5)')
        plt.scatter(jaccards,plot)
        plt.show()

    if arg == 'card':
        num_cards = 25
        cardinalities = np.logspace(1,5,num_cards)
        plot = np.zeros(num_cards)
        for i in range(num_cards):
            print('Starting '+str(i+1)+' out of '+str(num_cards))
            card = int(cardinalities[i])
            # do 5 trials for each cardinality
            trials = 10
            results = np.zeros(trials)
            for j in range(trials):
                h1 = HLL(12)
                h2 = HLL(12)
                a,b,exp_jaccard = Random_gen.generate_reads(0.1, card, card, 40)

                for s in a:
                    h1.insert(s)
                for s in b:
                    h2.insert(s)

                union = Similarity.union(h1,h2).cardinality()
                intersection = Similarity.intersection_inclusion_exclusion(h1,h2)

                obs_jaccard = intersection/union
                error = 100*(obs_jaccard - exp_jaccard)/(1+exp_jaccard)
                results[j] = error

            plot[i] = np.mean(results)

        print(plot)
        plt.xscale('log')
        plt.title('Jaccard (set at 0.01) Accuracy for reads of length 40')
        plt.xlabel('Cardinality')
        plt.ylabel('% Error (mean of 10)')
        plt.scatter(cardinalities,plot)
        plt.show()

if (__name__ == '__main__'):
    arg = sys.argv[1]
    main(arg)
