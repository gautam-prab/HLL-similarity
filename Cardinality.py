"""
Cardinality.py: Cardinality estimation implementation for our project
Functions in this file are implemented as class methods in HLL.py, so you shouldn't need to use this file.

Using cardinality estimation algorithms from:
Ertl, O. (2017). New cardinality estimation algorithms for HyperLogLog sketches. ArXiv.
"""
import numpy as np
import math

"""
Cardinality.getMultiplicity():
Get count vector from HLL registers. This is a helper method for the maximum likelihood estimation
Comes from Algorithm 4 in Ertl paper

Input: hll, a HLL object
Output: C, a list of integers, where C[i] represents the count of i in hll's registers.
"""
def getMultiplicity(hll):
    registers = hll.getRegisters()
    m = len(registers)

    one_position = hll.q
    C = []
    # Fills C with times K[i] appears in K
    for q in range(0, one_position + 1):
        C.append(0)

    # K[i] is a value in the range [1, q+1]
    for i in range(1,m):
        C[registers[i]] += 1
    return C

"""
Cardinality.estimateCardinality():
Calculate maximum likelihood estimator for cardinality of HLL
Comes from Algorithm 8 in Ertl paper

Input: counts, a count vector from Cardinality.getMultiplicity()
Output: a float representing the cardinality estimate for counts
"""
def estimateCardinality(counts):
    numCounts = len(counts) - 1
    q = numCounts - 1

    m = sum(counts)

    if counts[numCounts] == m:
        return math.inf

    #find the min K that has a nonzero count
    k_min = next((x for x, val in enumerate(counts) if val), None)
    k_prime_min = max(k_min, 1)
    k_max = next((x for x, val in enumerate(reversed(counts)) if val), None)
    k_max = numCounts - k_max
    k_prime_max = min(k_max, q)

    z = 0
    for k in range(k_prime_max, k_prime_min-1, -1):
        z = 0.5 * (z) + counts[k]
    z = z * 2 ** (-1 * k_prime_min)

    c = counts[q+1]
    if q >= 1:
        c = c + counts[k_prime_max]
    g_prev = 0
    a = z + counts[0]
    b = z + counts[q+1] * (2 ** (-1 * q))
    m_prime = m - counts[0]
    if b <= (1.5 * a):
        x = m_prime / (0.5 * b + a)  # weak lower bound
    else:
        x = (m_prime / (b)) * np.log(1 + b/a)  # strong lower bound

    delta_x = x
    delta = (10 ** -2) / np.sqrt(m)
    g_prev = 0
    while delta_x > x * delta: # maximum likelihood by Taylor approximation
        ka = math.floor(np.log2(x))
        x_prime = x * (2 ** (-1*max(k_prime_max+1, ka+1)))
        x_dprime = x_prime * x_prime
        #taylor approximation
        h = x_prime - x_dprime / 3 + (x_dprime * x_dprime) * (1/45 - x_dprime/472.5)
        for k in range(ka-1, k_prime_max-1, -1):
            print(x / (2 ** (k + 2)))
            top = x_prime + h*(1-h)
            bottom = x_prime + (1-h)
            h = top/bottom
            x_prime = 2 * x_prime
        g = c * h
        for k in range(k_prime_max - 1, k_prime_min-1, -1):
            top = x_prime + h * (1-h)
            bottom = x_prime + (1-h)
            h = top/bottom
            g = g + counts[k] * h
            x_prime = 2 * x_prime
        g = g + x * a
        if g > g_prev and m_prime >= g:
            delta_x = delta_x * (m_prime - g) / (g - g_prev)
        else:
            delta_x = 0
        x = x + delta_x
        g_prev = g
    return m * x

"""
Cardinality.simpleCardinality():
Calculate cardinality of HLL by bias-corrected harmonic mean
Comes from Algorithm 5 in Ertl paper

Input: registers, a list representing the registers of the HyperLogLog
Output: a float representing the cardinality estimate for registers
"""
def simpleCardinality(registers):
    m = len(registers)
    alpha_m = 0.7213/(1+1.079/m) # bias correction factor

    # raw estimation via harmonic mean
    n_raw = alpha_m * (m ** 2) / sum([(2 ** -x) for x in registers])

    # for sufficiently small cardinalities, we can just count number of 0 registers
    if (n_raw <= m*5/2): #small correction
        c0 = registers.count(0)
        if c0 != 0:
            return m*np.log2(m/c0)
        else:
            return n_raw
    elif (n_raw <= (2**64)/30):
        return n_raw
    else:
        # for sufficiently large cardinalities, we take into account fact that registers are saturated
        return -(2**64)*np.log2(1-n_raw/(2 ** 64))  #large correction
