from HLL import HLL
from Random_Generators import Rangen_jaccard

import numpy as np
import matplotlib.pyplot as plt
import sys

"""
PlotCardinalityAccuracy.main():
Calculates and plots cardinality percent error over 80 sets of different cardinalities over a logspace range
Inputs: num_cards - the number of sets
base_start - Smallest cardinality in range
base_stop - Largest cardinality in range
num_trials - number of trials cardinality calculations done for each set
"""

def main(num_cards, base_start, base_stop, num_trials, read_length):
    cardinalities = np.logspace(base_start, base_stop, num_cards)
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        results = np.zeros(num_trials)
        for j in range(num_trials):
            h = HLL(12)
            for k in range(card):
                # Generates a random sequence of length read_length and inserts it into HLL
                h.insert(Rangen_jaccard.generate_random_string(read_length))

            obs_cardinality = h.cardinality()  # Retrieves the cardinality

            results[j] = obs_cardinality

        plot[i] = 100 * (np.mean(results) - card) / card  # Plots percent error of calculated cardinality vs. expected value

    print('Mean Errors:')
    print(plot)
    plt.xscale('log')
    plt.title('Cardinality Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.ylim(-100,100)
    plt.scatter(cardinalities, plot)
    plt.show()

if (__name__ == '__main__'):
    # Hyperparameters used for this project: num_c = 80, card_s = 1, card_sp = 6.7, num_t = 10, read_length = 40
    num_c = int(sys.argv[1])
    card_s = float(sys.argv[2])
    card_sp = float(sys.argv[3])
    num_t = int(sys.argv[4])
    read_length = 40
    print(num_c, card_s, card_sp, num_t, read_length)
    main(num_c, card_s, card_sp, num_t, read_length)
