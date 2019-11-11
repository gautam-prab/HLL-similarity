import numpy as np
import math
'''Algorithm 4'''
def getMultiplicity(hll):
    registers = hll.getRegisters()
    m = len(registers)
    one_position = hll.q
    C = []
    # Fills C with times K[i] appears in K
    for q in range(0, one_position + 1):
        C.append(q)

    # K[i] is a value in the range [1, q+1]
    for i in range(1,m):
        if registers[i] not in C:
            C[registers[i]] = 1
        else:
            C[registers[i]] += 1
    return C


def estimateCardinality(counts, registers):
    numCounts = len(counts)
    q = numCounts - 1
    m = len(registers)

    if counts[numCounts] = m
        return math.inf

    #min K that has a nonzero count
    k_min = next(x for x, val in enumerate(counts) if val > 0)
    k_prime_min = max(k_min, 1)
    k_max = next(x for x, val in enumerate(reversed(counts) if val > 0)
    k_prime_max = min(k_max, q)

    z = 0
    for k in range(k_prime_min, k_prime_max):
        z = 0.5(z) + counts[k]
    z = z * 2^(-1 * k_prime_min)
    c = counts[numCounts]
    if q >= 1:
        c = c + couns[k_prime_max]
    g_prev = 0
    a = z + counts[0]
    b = z + counts[numCounts] * (2 ^ (-1 * q))
    m_prime = m - counts[0]
    if b <= (1.5 * a):
        x = m_prime / (0.5 * b + a)  # weak lower bound
    else:
        x = m_prime / (b * np.log(1 + b/a))  # strong lower bound

    delta_x = x
    delta = (10 ^ -2) / np.sqrt(m)
    while delta_x > x * delta:
        ka = 2 + floor(np.log2(x))
        x_prime = x * 2 ^(-max(k_max, ka)-1)
        x_dprime = x_prime * x_prime
        #taylor approximation
        h = x_prime - x_dprime/3 + (x_dprime * x_dprime) * (1/45 - x_dprime/472.5)
        for k in range(k_prime_max, ka -1):
            top = x_prime + h*(1-h)
            bottom = x_prime + (1-h)
            h = top/bottom
            x_prime = 2 * x_prime
        g = c * h
        for k in range(k_prime_min, k_prime_max - 1):
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
