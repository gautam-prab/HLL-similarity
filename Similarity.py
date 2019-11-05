import HLL, Cardinality
# Union and Intersection estimation implementation for our project
# Using cardinality estimation algorithms from:
    # Ertl, O. (2017). New cardinality estimation algorithms for HyperLogLog sketches. ArXiv.

def union(hll1, hll2):
    """Calculate the cardinality of the union of two sets"""
    # for now, assume sets have the same num. of registers; will fix later
    assert isinstance(hll1, HLL.HLL)
    assert isinstance(hll2, HLL.HLL)
    r1 = hll1.getRegisters()
    r2 = hll2.getRegisters()

    r_new = [0 for i in range(len(r1))]
    for i in range(len(r1)):
        r_new[i] = max(r1[i],r2[i])

    new_hll = HLL.HLL(hll1.p, r_new)
    return new_hll.cardinality()

def intersection(hll1, hll2):
    """Calculate the cardinality of the intersection of two sets"""
    assert isinstance(hll1, HLL.HLL)
    assert isinstance(hll2, HLL.HLL)
