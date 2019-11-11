import numpy as np
'''Algorithm 4'''
def getMultiplicity(hll):
    m = len(K)
    C = {}
    # Fills C with times K[i] appears in K
    # K[i] is a value in the range [1, q+1]
    for i in range(1, m):
        if K[i] not in C:
            C[K[i]] = 1
        else:
            C[K[i]] += 1
    return C


def estimateCardinality(counts):
    return

# Compute simple cardinality estimation, Algorithm 5 from Ertl
def simpleCardinality(registers):
    m = len(registers)
    alpha_m = 0.7213/(1+1.079/m)

    n_raw = alpha_m * (m ** 2) / sum([(2 ** -x) for x in registers])

    if (n_raw <= m*5/2): #small correction
        c0 = registers.count(0)
        if c0 != 0:
            return m*np.log2(m/c0)
        else:
            return n_raw
    elif (n_raw <= (2**64)/30):
        return n_raw
    else:
        return -(2**64)*np.log2(1-n_raw/(2 ** 64))  #large correction
