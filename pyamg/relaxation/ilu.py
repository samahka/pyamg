"""Incomplete LU factorization for sparse matrices"""

from warnings import warn

import numpy as np
import scipy as sp
from scipy import sparse

from pyamg import amg_core
from scipy.linalg import lapack as la

def make_system(A, formats=None):
     """
    Return A suitable for ILU or raise an exception
    Parameters
    ----------
    A : {sparse-matrix}
        n x n system
    formats: {'csr', 'csc', 'bsr', 'lil', 'dok',...}
        desired sparse matrix format
        default is no change to A's format
    Returns
    -------
    (A), where A is in the CSR sparse-matrix format.
    Notes
    -----
    Does some rudimentary error checking on the system,
    such as checking for compatible dimensions and checking
    for compatible type, i.e. float or complex.
    Examples
    --------
    >>> from pyamg.relaxation.relaxation import make_system
    >>> from pyamg.gallery import poisson
    >>> import numpy as np
    >>> A = poisson((10,10), format='csr')
    >>> (A) = make_system(A,formats=['csc'])
    """

    if formats is None:
        pass
    elif formats == ['csr']:
        if sparse.isspmatrix_csr(A):
            pass
        elif sparse.isspmatrix_bsr(A):
            A = A.tocsr()
        else:
            warn('implicit conversion to CSR', sparse.SparseEfficiencyWarning)
            A = sparse.csr_matrix(A)
    else:
        if sparse.isspmatrix(A) and A.format in formats:
            pass
        else:
            A = sparse.csr_matrix(A).asformat(formats[0])

    M, N = A.shape

    if M != N:
        raise ValueError('expected square matrix')

    return A

def ilu_k(A):

    





