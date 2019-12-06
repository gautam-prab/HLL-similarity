from HLL import HLL
import Similarity
from Random_Generators import Rangen_jaccard, Rangen_forbes, Rangen_sorenson_dice

import sys

import numpy as np
import matplotlib.pyplot as plt
import math

"""USAGE:
    python PlotSimilarityAccuracy.py jac to modulate jaccard
    python PlotSimilarityAccuracy.py card to modulate cardinality
"""
def intersection():
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
            # union = Similarity.union(h1,h2).cardinality()
            intersection = Similarity.intersection(h1,h2)
            #obs_jaccard = intersection/union
            num_overlapped = math.ceil(0.02 * (card * 2) / (0.02 + 1))
            error = 100*(intersection - num_overlapped)/(num_overlapped)
            results[j] = error
        plot[i] = np.mean(results)
    print(plot)
    plt.xscale('log')
    plt.title('Intersection Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.scatter(cardinalities,plot)
    plt.show()

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
            intersection = Similarity.intersection(h1,h2)
            obs_jaccard = intersection/union
            error = 100*(obs_jaccard - exp_jaccard)/(exp_jaccard)
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
            a,b,exp_forbes = Rangen_forbes.generate_reads(0.1, card, card, 40)
            for s in a:
                h1.insert(s)
            for s in b:
                h2.insert(s)
            union = Similarity.union(h1,h2).cardinality()
            intersection = Similarity.intersection(h1,h2)
            b_size, c_size, a_size = Similarity.getJointEstimators(h1, h2)
            obs_forbes = (intersection * union) / ((intersection * union) + 3 / 2 * (b_size * c_size))
            error = 100*(obs_forbes - exp_forbes)/(exp_forbes)
            results[j] = error
        plot[i] = np.mean(results)
        print(plot[i])

    print(plot)
    plt.xscale('log')
    plt.title('Forbes (set at 0.01) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.ylim(-100, 100)
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
            intersection = Similarity.intersection(h1,h2)
            obs_sd = 2 * intersection / (h1.cardinality + h2.cardinality)
            error = 100*(obs_sd - exp_sd)/(exp_sd)
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
    elif arg == "intersection":
        intersection()


if (__name__ == '__main__'):
    arg = sys.argv[1]
    main(arg)
