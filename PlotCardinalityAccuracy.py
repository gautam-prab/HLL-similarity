from HLL import HLL
import Random_gen

import numpy as np
import matplotlib.pyplot as plt

def main():

    num_cards = 50
    cardinalities = np.logspace(1,5,num_cards)
    plot = np.zeros(num_cards)
    for i in range(num_cards):
        print('Starting '+str(i+1)+' out of '+str(num_cards))
        card = int(cardinalities[i])
        # do 10 trials for each cardinality
        trials = 10
        results = np.zeros(trials)
        for j in range(trials):
            h = HLL(12)
            for k in range(card):
                h.insert(Random_gen.generate_random_string(40))

            obs_cardinality = h.cardinality()

            results[j] = obs_cardinality

        plot[i] = 100*(np.mean(results)-card)/card

    print(plot)
    plt.xscale('log')
    plt.title('Cardinality Accuracy for reads of length 40')
    plt.xlabel('Cardinality')
    plt.ylabel('% Error (mean of 10)')
    plt.scatter(cardinalities,plot)
    plt.show()

if (__name__ == '__main__'):
    main()
