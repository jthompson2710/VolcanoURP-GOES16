#----------------------------------------------------------------------------#
# Created on Mon Dec 14 13:12:33 2020
# Last edit September 14, 2024
# James Thompson - University of Texas at Austin, University of Pittsburgh
# For use by the University of West Indies
#----------------------------------------------------------------------------#
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

def baseline_als_optimized(y, lam, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z