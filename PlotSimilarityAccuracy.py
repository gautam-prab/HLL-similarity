from HLL import HLL
import Similarity
from Random_Generators import Rangen_jaccard, Rangen_forbes, Rangen_sorenson_dice

import sys

import numpy as np
import matplotlib.pyplot as plt
import math

"""USAGE:
    python PlotSimilarityAccuracy.py jac to modulate jaccard
    python PlotSimilarityAccuracy.py intersection to modulate intersection
    python PlotSimilarityAccuracy.py forbes to modulate forbes
    python PlotSimilarityAccuracy.py sd to modulate sorenson-dice
"""

"""
PlotSimilarityAccuracy.intersection():
Calculate and plot the percent error of intersection between two HLLs of various cardinalities
Based on an expected Jaccard value, using a simulated read generator based on a pre-set Jaccard
"""
def intersection():
    num_cards = 50
    cardinalities = np.logspace(2, 6.7, num_cards) # Creates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting ' + str(i + 1) + ' out of ' + str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        trials = 10
        results = np.zeros(trials)
        for j in range(trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Random read (length 40) generator based on the cardinalities and the expected Jaccard value
            a, b, exp_jaccard, forbes, sd = Rangen_jaccard.generate_reads(0.02, card, card, 40)
            for s in a:
                h1.insert(s)
            for s in b:
                h2.insert(s)
            intersection = Similarity.intersection(h1,h2) # Calculation of the intersection between 2 HLLs
            num_overlapped = math.ceil(0.02 * (card * 2) / (0.02 + 1)) # Calculation of the expected intersection based on the Jaccard formula
            error = 100 * (intersection - num_overlapped) / num_overlapped  # Percent error calculation
            results[j] = error
        plot[i] = np.mean(results)
    print(plot)  # Print out percent errors
    plt.xscale('log')
    plt.title('Intersection Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.ylim(-100, 100)
    plt.scatter(cardinalities, plot)
    plt.show()

"""
PlotSimilarityAccuracy.jaccard():
Calculate and plot the percent error of Jaccard coefficient between two HLLs of various cardinalities
Based on an expected Jaccard value, using a simulated read generator based on a pre-set Jaccard
"""
def jaccard():
    num_cards = 50
    cardinalities = np.logspace(2, 6.7, num_cards)  # Generates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        trials = 10
        results = np.zeros(trials)
        for j in range(trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Generate reads of length 40 based on the expected Jaccard, 0.02
            a, b, exp_jaccard, forbes, sd = Rangen_jaccard.generate_reads(0.02, card, card, 40)
            for s in a:
                h1.insert(s)
            for s in b:
                h2.insert(s)
            # Union and intersection calculations for Jaccard calculations
            union = Similarity.union(h1,h2).cardinality()
            intersection = Similarity.intersection(h1,h2)
            obs_jaccard = intersection/union
            error = 100 * (obs_jaccard - exp_jaccard) / exp_jaccard # Expected Jaccard is 0.02, error calculation
            results[j] = error
        plot[i] = np.mean(results)
    print(plot) # Percent Errors being displayed
    plt.xscale('log')
    plt.title('Jaccard (set at 0.02) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.ylim(-100, 100)
    plt.scatter(cardinalities,plot)
    plt.show()

"""
PlotSimilarityAccuracy.forbes():
Calculate and plot the percent error of Forbes coefficient between two HLLs of various cardinalities
Based on an expected Forbes value, using a simulated read generator based on a pre-set Forbes
"""
def forbes():
    num_cards = 50
    cardinalities = np.logspace(2, 6.7, num_cards) # Generates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        trials = 10
        results = np.zeros(trials)
        for j in range(trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Generate reads of length 40 based on the expected Forbes, 0.1
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

    print(plot)
    plt.xscale('log')
    plt.title('Forbes (set at 0.1) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.ylim(-100, 100)
    plt.scatter(cardinalities,plot)
    plt.show()

"""
PlotSimilarityAccuracy.sd():
Calculate and plot the percent error of Sorensen-Dice coefficient between two HLLs of various cardinalities
Based on an expected Sorensen-Dice value, using a simulated read generator based on a pre-set Sorensen-Dice
"""
def sd():
    num_cards = 50
    cardinalities = np.logspace(2, 6.7, num_cards) # Generates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        trials = 10
        results = np.zeros(trials)
        for j in range(trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Generate reads of length 40 based on the expected Forbes, 0.04
            a,b,exp_sd = Rangen_sorenson_dice.generate_reads(0.04, card, card, 40)
            for s in a:
                h1.insert(s)
            for s in b:
                h2.insert(s)
            # Calculation of intersection between two HLLs, necessary for SD calculation
            intersection = Similarity.intersection(h1,h2)
            obs_sd = 2 * intersection / (h1.cardinality() + h2.cardinality()) # Calculation of expected Sorensen-Dice
            error = 100 * (obs_sd - exp_sd) / exp_sd # Percent error calculation for SD
            results[j] = error
        plot[i] = np.mean(results)
    print(plot) # Print out all percent errors
    plt.xscale('log')
    plt.title('Sorenson-Dice (set at 0.04) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.scatter(cardinalities,plot)
    plt.show()

"""
main()
Calls the appropriate function based on the command line argument received
Input: The command-line argument
"""
def main(arg):
    # Function callers based on command line argument
    if arg == "jaccard":
        jaccard()
    elif arg == "forbes":
        forbes()
    elif arg == "sd":
        sd()
    elif arg == "intersection":
        intersection()


if (__name__ == '__main__'):
    arg = sys.argv[1]
    main(arg)
