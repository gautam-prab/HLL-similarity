import Cardinality
# HyperLogLog implementation for our project
# Using cardinality estimation algorithms from:
    # Ertl, O. (2017). New cardinality estimation algorithms for HyperLogLog sketches. ArXiv.

class HLL:
    """HyperLogLog implementation"""

    def __init__(self, p):
        """Initialize HLL with 2^p registers"""
        self.p = p
        self.q = 64 - p # assuming 64-bit hash value
        self.m = 2 ** p

    def insert(self, read):
        """Insert read into the HLL"""
        return

    def cardinality(self):
        C = Cardinality.getMultiplicity()
        Cardinality.estimateCardinality(C)
        return
