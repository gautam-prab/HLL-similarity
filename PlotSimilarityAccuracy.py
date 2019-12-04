from HLL import HLL
import Similarity
from Random_Generators import Rangen_jaccard, Rangen_forbes, Rangen_sorenson_dice

import sys

import numpy as np
import matplotlib.pyplot as plt

"""USAGE:
    python PlotSimilarityAccuracy.py jac to modulate jaccard
    python PlotSimilarityAccuracy.py card to modulate cardinality
"""

def jaccard():
    num_cards = 50
    cardinalities = np.logspace(2, 6.7, num_cards)
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
            a,b,exp_jaccard,forbes,sd = Rangen_jaccard.generate_reads(0.02, card, card, 40)
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

def forbes():
    num_cards = 50
    cardinalities = np.logspace(2, 6.7, num_cards)
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
            a,b,exp_forbes = Rangen_forbes.generate_reads(0.0008, card, card, 40)
            for s in a:
                h1.insert(s)
            for s in b:
                h2.insert(s)
            union = Similarity.union(h1,h2).cardinality()
            intersection = Similarity.intersection_inclusion_exclusion(h1,h2)
            obs_forbes = intersection/union
            error = 100*(obs_forbes - exp_forbes)/(1+exp_forbes)
            results[j] = error
        plot[i] = np.mean(results)
    print(plot)
    plt.xscale('log')
    plt.title('Forbes (set at 0.01) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.scatter(cardinalities,plot)
    plt.show()
    
def sd():
    num_cards = 50
    cardinalities = np.logspace(2, 6.7, num_cards)
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
            a,b,exp_sd = Rangen_sorenson_dice.generate_reads(0.04, card, card, 40)
            for s in a:
                h1.insert(s)
            for s in b:
                h2.insert(s)
            union = Similarity.union(h1,h2).cardinality()
            intersection = Similarity.intersection_inclusion_exclusion(h1,h2)
            obs_sd = intersection/union
            error = 100*(obs_sd - exp_sd)/(1+exp_sd)
            results[j] = error
        plot[i] = np.mean(results)
    print(plot)
    plt.xscale('log')
    plt.title('Sorenson-Dice (set at 0.01) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.scatter(cardinalities,plot)
    plt.show()

def main(arg):
    if arg == "jaccard":
        jaccard()
    elif arg == "forbes":
        forbes()
    elif arg == "sorenson-dice":
        sd()

if (__name__ == '__main__'):
    arg = sys.argv[1]
    main(arg)
