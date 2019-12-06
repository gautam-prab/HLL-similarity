# HLL-similarity
Set similarity methods for HyperLogLog sketches for genomic comparisons

JHU computational genomics project 2019

1. To produce the cardinality graph, run the following line in the Terminal:

        python3 PlotCardinalityAccuracy.py

In this particular files, there are the following hyperparameters:
* num_cards - the number of cardinalities that will be used
* trials - the number of trials that we will average across
The above code takes around 10 hrs, for the hyperparameters we used.



2. To produce each of the similarity graphs, run the following line in the Terminal:

        python3 PlotCardinalityAccuracy.py "similarity-metric"

where you replace "similarity-metric" with jaccard/forbes/sorenson-dice/intersection

In this particular files, there are the following hyperparameters:
* num_cards - the number of cardinalities that will be used
* trials - the number of trials that we will average across
* whatever the value the similarity metric has
The above code takes around 13 hrs for each, for the hyperparameters we used.
