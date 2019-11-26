import Cardinality
from HLL import HLL
import numpy as np
import math

from Random_Generators import Rangen_jaccard
from scipy.optimize import fmin_bfgs
from scipy.optimize import minimize
# Union and Intersection estimation implementation for our project
# Using cardinality estimation algorithms from:
    # Ertl, O. (2017). New cardinality estimation algorithms for HyperLogLog sketches. ArXiv.

def union(hll1, hll2):
    """Calculate the HLL corresponding to the union of two sets"""
    # for now, assume hlls have the same num. of registers; will fix later
    assert isinstance(hll1, HLL)
    assert isinstance(hll2, HLL)
    r1 = hll1.getRegisters()
    r2 = hll2.getRegisters()

    r_new = [0 for i in range(len(r1))]
    for i in range(len(r1)):
        r_new[i] = max(r1[i],r2[i])

    new_hll = HLL(hll1.p, r_new)
    return new_hll

def intersection_inclusion_exclusion(hll1, hll2):
    """Calculate cardinality of the intersection through the inclusion-exclusion principle"""
    return hll1.cardinality() + hll2.cardinality() - union(hll1,hll2).cardinality()

def main():
    h1 = HLL(8)
    h2 = HLL(8)

    # generate 100,000 reads for each
    # make them have a Jaccard similarity of 0.05 (1/20 reads are shared)
    a,b,jac = Rangen_jaccard.generate_reads(0.05, 100000, 100000, 50)

    for r in a:
        h1.insert(r)
    for r in b:
        h2.insert(r)

    print('H1 Cardinality: '+str(h1.cardinality()))
    print('H2 Cardinality: '+str(h2.cardinality()))

    print('Expected jaccard: '+str(jac))
    u = union(h1,h2).cardinality()
    i = intersection_inclusion_exclusion(h1,h2)
    print('Union: '+str(u))
    print('Intersect: '+str(i))
    print('Received Jaccard: '+str(i/u))

def intersection(hll1, hll2):
    """Calculate the cardinality of the intersection of two sets"""
    # for now, assume hlls have same num. of registers
    assert isinstance(hll1, HLL)
    assert isinstance(hll2, HLL)
    r1 = hll1.getRegisters()
    r2 = hll2.getRegisters()

    q = hll1.q
    m = hll1.m

    C1_less = np.zeros(q+2)
    C2_less = np.zeros(q+2)
    C1_greater = np.zeros(q+2)
    C2_greater = np.zeros(q+2)
    C_equal = np.zeros(q+2)

    for i in range(len(r1)):
        if r1[i] < r2[i]:
            C1_less[r1[i]] += 1
            C2_greater[r2[i]] += 1
        elif r1[i] > r2[i]:
            C1_greater[r1[i]] += 1
            C2_less[r2[i]] += 1
        else:
            C_equal[r1[i]] += 1

    lambda_ax = Cardinality.estimateCardinality(C1_less + C_equal + C1_greater)
    lambda_bx = Cardinality.estimateCardinality(C2_greater + C_equal + C2_less)

    if C1_less[0] + C_equal[0] + C2_less[0] == m:
        # at least one HLL is 0 at every position
        return 0

    lambda_abx = Cardinality.estimateCardinality(C1_greater + C_equal + C2_greater)

    phi_a = np.log(max(1, lambda_abx-lambda_bx))
    phi_b = np.log(max(1, lambda_abx-lambda_ax))
    phi_x = np.log(max(1, lambda_ax+lambda_bx-lambda_abx))

    optimize_condition = True
    delta = 0.01/np.sqrt(m)

    args = (C1_less, C2_less, C_equal, C1_greater, C2_greater, q, m)
    # results = fmin_bfgs(calculate_log_likelihood, np.array([phi_a+50, phi_b+50, phi_x+50]), calculate_ll_gradient, args=args, gtol=delta, disp=True, epsilon = 1.5e-2)
    results = minimize(calculate_log_likelihood, np.array([phi_a+1, phi_b+1, phi_x+1]), method='L-BFGS-B', args=args, options={'gtol':delta, 'disp':False})

    phi_a = results.x[0] # |A \ B|
    phi_b = results.x[1] # |B \ A|
    phi_x = results.x[2] # |A n B|

    return np.exp(phi_x)

# Equation 72 from Ertl paper
def calculate_log_likelihood(phi, C1_less, C2_less, C_equal, C1_greater, C2_greater, q, m):
    phi_a = phi[0]
    phi_b = phi[1]
    phi_x = phi[2]

    f = 0
    (xa,ya,za) = calculate_xyz(phi_a, q, m)
    (xb,yb,zb) = calculate_xyz(phi_b, q, m)
    (xx,yx,zx) = calculate_xyz(phi_x, q, m)
    for k in range(1,q+1):
        (xak,yak,zak) = (xa[k],ya[k],za[k])
        (xbk,ybk,zbk) = (xb[k],yb[k],zb[k])
        (xxk,yxk,zxk) = (xx[k],yx[k],zx[k])
        f += C1_less[k]*np.log(zxk+yxk*zak) + C2_less[k]*np.log(zxk+yxk*zbk) + C1_greater[k]*np.log(zak) + C2_greater[k]*np.log(zbk) + C_equal[k]*np.log(zxk+yxk*zak*zbk)
    f += C1_greater[q+1]*np.log(za[q]) + C2_greater[q+1]*np.log(zb[q]) + C_equal[q+1]*np.log(zx[q]+yx[q]*za[q]*zb[q])
    for k in range(0,q+1):
        f -= (C1_less[k] + C_equal[k] + C1_greater[k])*xa[k] + (C2_less[k] + C_equal[k] + C2_greater[k])*xb[k] + (C1_less[k] + C_equal[k] + C2_less[k])*xx[k]
    return -f

# Not using LL gradient – there is an error in it and it's more accurate to just calculate gradient by finite differences

# def calculate_ll_gradient(phi, C1_less, C2_less, C_equal, C1_greater, C2_greater, q, m):
#     phi_a = phi[0]
#     phi_b = phi[1]
#     phi_x = phi[2]
#
#     da = 0
#     db = 0
#     dx = 0
#     (xa,ya,za) = calculate_xyz(phi_a, q, m)
#     (xb,yb,zb) = calculate_xyz(phi_b, q, m)
#     (xx,yx,zx) = calculate_xyz(phi_x, q, m)
#     for k in range(1,q+1):
#         da += C1_less[k]*(yx[k]*xa[k]*ya[k])/(zx[k]+yx[k]*za[k]) + C1_greater[k]*(xa[k]*ya[k])/(za[k]) + C_equal[k]*(yx[k]*xa[k]*ya[k]*zb[k])/(zx[k]+yx[k]*za[k]*zb[k])
#         db += C2_less[k]*(yx[k]*xb[k]*yb[k])/(zx[k]+yx[k]*zb[k]) + C2_greater[k]*(xb[k]*yb[k])/(zb[k]) + C_equal[k]*(yx[k]*za[k]*xb[k]*yb[k])/(zx[k]+yx[k]*za[k]*zb[k])
#         dx += C1_less[k]*(xx[k]*yx[k]*ya[k])/(zx[k]+yx[k]*za[k]) + C2_less[k]*(xx[k]*yx[k]*yb[k])/(zx[k]+yx[k]*zb[k]) + C_equal[k]*(xx[k]*yx[k]*(ya[k]+za[k]*yb[k]))/(zx[k]+yx[k]*za[k]*zb[k])
#     da += C1_greater[q+1]*(xa[q]*ya[q])/(za[q]) + C_equal[q+1]*(yx[q]*xa[q]*ya[q]*zb[q])/(zx[q]+yx[q]*za[q]*zb[q])
#     db += C2_greater[q+1]*(xb[q]*yb[q])/(zb[q]) + C_equal[q+1]*(yx[q]*za[q]*xb[q]*yb[q])/(zx[q]+yx[q]*za[q]*zb[q])
#     dx += C_equal[q+1]*(xx[q]*yx[q]*(ya[q]+za[q]*yb[q]))/(zx[q]+yx[q]*za[q]*zb[q])
#     for k in range(0,q+1):
#         da -= (C1_less[k] + C_equal[k] + C1_greater[k])*xa[k]
#         db -= (C2_less[k] + C_equal[k] + C2_greater[k])*xb[k]
#         dx -= (C1_less[k] + C_equal[k] + C2_less[k])*xx[k]
#     return np.array([da, db, dx])

def calculate_xyz(phi, q, m):
    x = np.zeros(q+1)
    z = np.zeros(q+1)
    y = np.zeros(q+1)
    for k in range(0,q+1):
        x[k] = np.exp(phi) / (m*(2 ** k))
        if x[k] < np.log(2):
            z[k] = -np.expm1(-x[k])
            y[k] = 1 - z[k]
        else:
            y[k] = np.exp(-x[k])
            z[k] = 1 - y[k]
    return (x,y,z)

def intersection_testing():
    h1 = HLL(8)
    h2 = HLL(8)

    jac = 0.1
    card_a = 100000
    card_b = 100000

    a,b,jac,_ = Rangen_jaccard.generate_reads(jac, card_a, card_b, 50)

    for r in a:
        h1.insert(r)
    for r in b:
        h2.insert(r)

    print('H1 Cardinality: '+str(h1.cardinality()))
    print('H2 Cardinality: '+str(h2.cardinality()))

    expected_i = math.ceil(jac * (card_a + card_b) / (jac + 1))

    print('Expected intersection: '+str(expected_i))
    i1 = intersection_inclusion_exclusion(h1,h2)
    print('Intersect IE: '+str(i1))
    i2 = intersection(h1,h2)
    print('Intersect MLE: '+str(i2))

if __name__ == "__main__":
    main()
    # intersection_testing()
