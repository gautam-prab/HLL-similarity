
'''Algorithm 4'''
def getMultiplicity(K):
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


'''Algorithm 8'''
def estimateCardinality(C):
    m = len(C)
    if C[q+1] = m
