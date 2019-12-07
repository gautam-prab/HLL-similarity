from HLL import HLL
import Similarity
from Random_Generators import Rangen_jaccard, Rangen_forbes, Rangen_sorensen_dice

import sys

import numpy as np
import matplotlib.pyplot as plt
import math

"""USAGE:
    python PlotSimilarityAccuracy.py jac to modulate jaccard
    python PlotSimilarityAccuracy.py intersection to modulate intersection
    python PlotSimilarityAccuracy.py forbes to modulate forbes
    python PlotSimilarityAccuracy.py sd to modulate sorensen-dice
    Accompanied with appropriate command line argument parameters, described at bottom
"""

"""
PlotSimilarityAccuracy.intersection():
Calculate and plot the percent error of intersection between two HLLs of various cardinalities
Based on an expected Jaccard value, using a simulated read generator based on a pre-set Jaccard
Input: num_cards - number of cardinalities to generate
base_start - Smallest cardinality in range
base_stop - Largest cardinality in range
num_trials - number of trials (intersection calculations) done for each set's cardinality
"""

# num_cards, range of cardinalities, number of trials
def intersection(num_cards, base_start, base_stop, num_trials, exp_jaccard, read_lengths):
    cardinalities = np.logspace(base_start, base_stop, num_cards) # Creates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting ' + str(i + 1) + ' out of ' + str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        results = np.zeros(num_trials)
        for j in range(num_trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Random read (length 40) generator based on the cardinalities and the expected Jaccard value
            a, b, exp_jaccard, forbes, sd = Rangen_jaccard.generate_reads(exp_jaccard, card, card, read_lengths)
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
Input: num_cards - number of cardinalities to generate
base_start - Smallest cardinality in range
base_stop - Largest cardinality in range
num_trials - number of trials (Jaccard calculations) done for each set's cardinality
"""
def jaccard(num_cards, base_start, base_stop, num_trials, exp_jaccard, read_lengths):
    cardinalities = np.logspace(base_start, base_stop, num_cards)  # Generates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        results = np.zeros(num_trials)
        for j in range(num_trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Generate reads of length 40 based on the expected Jaccard, 0.02
            a, b, exp_jaccard, forbes, sd = Rangen_jaccard.generate_reads(exp_jaccard, card, card, read_lengths)
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
    plt.scatter(cardinalities, plot)
    plt.show()

"""
PlotSimilarityAccuracy.forbes():
Calculate and plot the percent error of Forbes coefficient between two HLLs of various cardinalities
Based on an expected Forbes value, using a simulated read generator based on a pre-set Forbes
Input: num_cards - number of cardinalities to generate
base_start - Smallest cardinality in range
base_stop - Largest cardinality in range
num_trials - number of trials (Forbes calculations) done for each set's cardinality
"""
def forbes(num_cards, base_start, base_stop, num_trials, exp_forbes, read_lengths):
    # num_cards = 50
    cardinalities = np.logspace(base_start, base_stop, num_cards) # Generates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        results = np.zeros(num_trials)
        for j in range(num_trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Generate reads of length 40 based on the expected Forbes, 0.1
            a,b,exp_forbes = Rangen_forbes.generate_reads(exp_forbes, card, card, read_lengths)
            for s in a:
                h1.insert(s)
            for s in b:
                h2.insert(s)
            union = Similarity.union(h1,h2).cardinality()
            intersection = Similarity.intersection(h1,h2)
            b_size, c_size, a_size = Similarity.getJointEstimators(h1, h2)
            obs_forbes = (intersection * union) / ((intersection * union) + 3 / 2 * (b_size * c_size))
            error = 100 * (obs_forbes - exp_forbes) / exp_forbes
            results[j] = error
        plot[i] = np.mean(results)

    print(plot)
    plt.xscale('log')
    plt.title('Forbes (set at 0.1) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.ylim(-100, 100)
    plt.scatter(cardinalities, plot)
    plt.show()

"""
PlotSimilarityAccuracy.sd():
Calculate and plot the percent error of Sorensen-Dice coefficient between two HLLs of various cardinalities
Based on an expected Sorensen-Dice value, using a simulated read generator based on a pre-set Sorensen-Dice
Input: num_cards - number of cardinalities to generate
base_start - Smallest cardinality in range
base_stop - Largest cardinality in range
num_trials - number of trials (Sorensen-Dice calculations) done for each set's cardinality
"""
def sd(num_cards, base_start, base_stop, num_trials, exp_sd, read_lengths):
    cardinalities = np.logspace(base_start, base_stop, num_cards) # Generates a range of cardinalities for set generation
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        results = np.zeros(num_trials)
        for j in range(num_trials):
            h1 = HLL(12)
            h2 = HLL(12)
            # Generate reads of length 40 based on the expected Sorensen-Dice, 0.04
            a,b,exp_sd = Rangen_sorensen_dice.generate_reads(exp_sd, card, card, read_lengths)
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
    print(plot)  # Print out all percent errors
    plt.xscale('log')
    plt.title('Sorensen-Dice (set at 0.04) Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.scatter(cardinalities, plot)
    plt.show()

"""
main()
Calls the appropriate function based on the command line argument received
Input: The command-line arguments
"""
def main(arg, num_c, card_s, card_sp, num_t, exp_const, read_lengths):
    # Function callers based on command line argument
    if arg == "jaccard":
        jaccard(num_c, card_s, card_sp, num_t, exp_const, read_lengths)
    elif arg == "forbes":
        forbes(num_c, card_s, card_sp, num_t, exp_const, read_lengths)
    elif arg == "sd":
        sd(num_c, card_s, card_sp, num_t, exp_const, read_lengths)
    elif arg == "intersection":
        intersection(num_c, card_s, card_sp, num_t, exp_const, read_lengths)


if (__name__ == '__main__'):
    # Hyperparameters we used for this project: num_c = 50, card_s = 2, card_sp = 6.7, num_t = 10
    # Specific hyperparameters we used for the different methods:
    # Intersection and Jaccard: exp_const = 0.02, read_lengths = 40
    # Forbes: exp_const = 0.1, read_lengths = 40
    # Sorensen-Dice: exp_const = 0.04, read_lengths = 40
    arg = sys.argv[1]
    num_c = int(sys.argv[2])
    card_s = float(sys.argv[3])
    card_sp = float(sys.argv[4])
    num_t = int(sys.argv[5])
    exp_const = float(sys.argv[6])
    read_lengths = 40
    print(arg, num_c, card_s, num_t, exp_const, read_lengths)
    main(arg, num_c, card_s, card_sp, num_t, exp_const, read_lengths)
