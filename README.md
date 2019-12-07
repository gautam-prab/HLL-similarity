# GNome: Sketching-Based Set Similarity for Genomic Comparisons

We tested different set similarity metrics and
1. Measure their ability to be accurately estimated by the HyperLogLog data structure
2. Measure their accuracy in estimating the Average Nucleotide Identity of a set of genomes

This is our final project for JHU's Computational Genomics class in 2019 by Justin Greene, Gautam Prabhu, Jonathan Wang, and David Yang.

### Creating a HyperLogLog
We've implemented the HyperLogLog data structure in `HLL.py`. Basic estimation algorithms can be found in `Cardinality.py` and `Similarity.py`.

1. As a basic example, the `main` function of `HLL.py` will create a random set of n reads, insert them into an HLL, and output its cardinality estimates. For 10,000 reads, use the following command:

        python3 HLL.py 10000

2. As an example of the estimators for union and intersection, the `main` function of `Similarity.py` will create two random sets of reads with a given Jaccard value. Both will be inserted into an HLL and all estimates will be outputted. For an example with 10,000 reads in both sets and a Jaccard of 0.04, use the following command:

        python3 Similarity.py 10000 10000 0.04

### Benchmarking HyperLogLog Accuracy

In order to show that our implementations produce reasonable accuracy, we generated plots of estimation error at a variety of cardinalities.

1. To see a graph of cardinality accuracy, we use  `PlotCardinalityAccuracy.py`. For an example with 10 cardinalities from 10<sup>1</sup> to 10<sup>5</sup>, averaging 2 trials per datapoint, run:

        python3 PlotCardinalityAccuracy.py 10 1 5 2

In our paper we used the parameters 80, 1, 6.7, 10. Running with these parameters took us around 10 hours.

2. To see a graph of similarity accuracy, we use  `PlotSimilarityAccuracy.py`. We can calculate accuracy for intersection, Jaccard coefficient, Sørensen-Dice coefficient, and Alroy's modified Forbes coefficient. For these accuracies with 10 cardinalities from 10<sup>1</sup> to 10<sup>4</sup>, averaging 2 trials per datapoint, run:

        python3 PlotSimilarityAccuracy.py [measure] 10 1 4 2 [expected-measure]

where [measure] is either from the set jaccard/forbes/sd/intersection and [expected-measure] is an expected value of the measure from 0-1 (we recommend 0.02 for Jaccard, 0.1 for Forbes, and 0.04 for Sørensen-Dice).

For instance, for the Jaccard function, run:

        python3 PlotSimilarityAccuracy.py jaccard 10 1 4 2 0.02

For our project, we ran all four of these with the following parameters:
* jaccard 50 2 6.7 10 0.02
* sd 50 2 6.7 10 0.04
* forbes 50 2 6.7 10 0.1
* intersection 50 2 6.7 10 0.02

Again, these took up to 13 hours on our machine.

### Running GNome

To generate genomic comparisons against a ground truth, we use `GNome.py`. To test this on a small dataset synthetic Illumina sequencing reads from 10 bacterial genomes at 0.5x coverage, use the command:

    python3 GNome.py Data/0.5x/

This will output similarity scores between +1 and -1. A score of +1 would mean that that similarity metric ranks genome similarity equivalently to the ground truth.

This command will take a few minutes to run (our system is not optimized for high efficiency). When prompted, you can save the HLL sketches for further analysis.

### References

In our project, we utilized some outside sources for the data found in this repository.

* We used the 64-bit mixing integer hash function from Thomas Wang's [Integer Hash Functions](https://gist.github.com/badboy/6267743)
* Synthetic Illumina reads at different coverages were generated from Manuel Holtgrewe's [Mason](https://www.seqan.de/apps/mason/)
* The ground truth for our dataset was generated using the Kostas Lab's [ANI Matrix](http://enve-omics.ce.gatech.edu/g-matrix/index)
