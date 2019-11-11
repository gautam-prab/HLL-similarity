import Cardinality
from HLL import HLL
import numpy as np

import Random_gen
# from scipy.optimize import minimize
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
    a,b,jac = Random_gen.generate_reads(0.05, 100000, 100000, 50)

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

if __name__ == "__main__":
    main()

# def intersection(hll1, hll2):
#     """Calculate the cardinality of the intersection of two sets"""
#     # for now, assume hlls have same num. of registers
#     assert isinstance(hll1, HLL.HLL)
#     assert isinstance(hll2, HLL.HLL)
#     r1 = hll1.getRegisters()
#     r2 = hll2.getRegisters()
#
#     q = hll1.q
#     m = hll1.m
#
#     C1_less = np.zeros(q+1)
#     C2_less = np.zeros(q+1)
#     C1_greater = np.zeros(q+1)
#     C2_greater = np.zeros(q+1)
#     C_equal = np.zeros(q+1)
#
#     for i in range(len(r1)):
#         if r1[i] < r2[i]:
#             C1_less[r1[i]] += 1
#             C2_greater[r2[i]] += 1
#         elif r1[i] > r2[i]:
#             C1_greater[r1[i]] += 1
#             C2_less[r2[i]] += 1
#         else:
#             c_equal += 1
#
#     lambda_ax = Cardinality.estimateCardinality(C1_less + C_equal + C1_greater)
#     lambda_bx = Cardinality.estimateCardinality(C2_greater + C_equal + C2_less)
#
#     if C1_less[0] + c_equal[0] + C2_less[0] == m:
#         return (lambda_ax, lambda_bx, 0)
#
#     lambda_abx = Cardinality.estimateCardinality(C1_greater + C_equal + C2_greater)
#
#     phi_a = log(max(1, lambda_abx-lambda_bx))
#     phi_b = log(max(1, lambda_abx-lambda_ax))
#     phi_x = log(max(1, lambda_ax+lambda_bx-lambda_abx))
#
#     optimize_condition = True
#     delta = 0.01/sqrt(m)
#
#     fun = lambda (phia, phib, phic): -1*calculate_log_likelihood(C1_less, C2_less, C1_greater, C2_greater, phia, phib, phix, q, m)
#     grad = lambda (phia, phib, phic): -1*calculate_ll_gradient(C1_less, C2_less, C1_greater, C2_greater, phia, phib, phix, q, m)
#
#     (phi_a, phi_b, phi_x) = minimize(fun, [phi_a, phi_b, phi_x], method='BFGS', jac=grad, tol=delta)
#
#     return (np.exp(phi_a), np.exp(phi_b), np.exp(phi_x))
#
# def calculate_log_likelihood(C1_less, C2_less, C_equal, C1_greater, C2_greater, phi_a, phi_b, phi_x, q, m):
#     f = 0
#     (xa,ya,za) = calculate_xyz(phi_a, m)
#     (xb,yb,zb) = calculate_xyz(phi_b, m)
#     (xx,yx,zx) = calculate_xyz(phi_x, m)
#     for k in range(1,q+1):
#         (xak,yak,zak) = (xa[k],ya[k],za[k])
#         (xbk,ybk,zbk) = (xb[k],yb[k],zb[k])
#         (xxk,yxk,zxk) = (xx[k],yx[k],zx[k])
#         f += C1_less[k]*log(zxk+yxk*zak) + C2_less[k]*log(zxk+yxk*zbk) + C1_greater[k]*log(zak) + C2_greater[k]*log(zbk) + C_equal[k]*log(zxk+yxk*zak*zbk)
#     f += C1_greater[q+1]*log(za[q]) + C2_greater[q+1]*log(zb[q]) + C_equal[q+1]*log(zx[q]+yx[q]*za[q]*zb[q])
#     for k in range(0,q+1):
#         (xak,yak,zak) = (xa[k],ya[k],za[k])
#         (xbk,ybk,zbk) = (xb[k],yb[k],zb[k])
#         (xxk,yxk,zxk) = (xx[k],yx[k],zx[k])
#         f -= (C1_less[k] + C_equal[k] + C1_greater[k])*xak + (C2_less[k] + C_equal[k] + C2_greater[k])*xbk + (C1_less[k] + C_equal[k] + C2_less[k])*xxk
#     return f
#
# def calculate_ll_gradient(C1_less, C2_less, C_equal, C1_greater, C2_greater, phi_a, phi_b, phi_x, q, m):
#     da = 0
#     db = 0
#     dx = 0
#     for k in range(1,q+1):
#         (xa,ya,za) = calculate_xyz(phi_a, k, m)
#         (xb,yb,zb) = calculate_xyz(phi_a, k, m)
#         (xx,yx,zx) = calculate_xyz(phi_a, k, m)
#         da += C1_less[k]*(yx*xa*ya)/(zx+yx*za) + C1_greater[k]*(xa*ya)/(za) + C_equal[k]*(yx*xa*ya*zb)/(zx+yx*za*zb)
#         db += C2_less[k]*(yx*xb*yb)/(zx+yx*zb) + C2_greater[k]*(xb*yb)/(zb) + C_equal[k]*(yx*za*xb*yb)/(zx+yx*za*zb)
#         dx += C1_less[k]*(xx*yx*ya)/(zx+yx*za) + C2_less[k]*(xx*yx*yb)/(zx+yx*zb) + C1_greater[k]*log(za) + C2_greater[k]*log(zb) + C_equal[k]*(xx*yx*(ya+za*yb))/(zx+yx*za*zb)
#     (xa,ya,za) = calculate_xyz(phi_a, q, m)
#     (xb,yb,zb) = calculate_xyz(phi_a, q, m)
#     (xx,yx,zx) = calculate_xyz(phi_a, q, m)
#     da += C1_greater[q+1]*(xa*ya)/(za) + C_equal[q+1]*(yx*xa*ya*zb)/(zx+yx*za*zb)
#     db += C2_greater[q+1]*(xb*yb)/(zb) + C_equal[q+1]*(yx*za*xb*yb)/(zx+yx*za*zb)
#     dx += C_equal[q+1]*(xx*yx*(ya+za*yb))/(zx+yx*za*zb)
#     for k in range(0,q+1):
#         (xa,ya,za) = calculate_xyz(phi_a, k, m)
#         (xb,yb,zb) = calculate_xyz(phi_a, k, m)
#         (xx,yx,zx) = calculate_xyz(phi_a, k, m)
#         da -= (C1_less[k] + C_equal[k] + C1_greater[k])*xa
#         db -= (C2_less[k] + C_equal[k] + C2_greater[k])*xb
#         dx -= (C1_less[k] + C_equal[k] + C2_less[k])*xx
#     return (da, db, dx)
#
# def calculate_xyz(phi, q, m):
#     x = np.zeros(q+1)
#     z = np.zeros(q+1)
#     y = np.zeros(q+1)
#     for k in range(0,q+1)
#         x[k] = np.exp(phi) / (m*(2 ** k))
#         z[k] = -np.expm1(-x[k])
#         y[k] = 1 - z[k]
#     return (x,y,z)
