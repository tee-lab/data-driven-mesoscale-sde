"""
Utilities for fitting vector and matrix polynomial functions of a given vector m.

VECTOR POLYNOMIALS
Given vectors m, a: R^d, the goal is to fit a function f: R^d -> R^d such that a = f(m).
Assuming that m and f(m) are rotationally symmetric, f should depend only on m as a whole,
and not on the individual components (mx, my).

We try to represent f(m) as a sparse polynomial of m, i.e, containing terms of the form:
m, |m| m, |m|^2 m, |m|^3 m, ...

VectorPoly is class that can represent, print and compute with vector polynomials of the above form.
VectorPolyFit implements sparse regression to fit a vector polynomial from given data.

MATRIX POLYNOMIALS
Given vector m: R^d and matrix B: R^(d*d), the goal is to fit a function G: R^d -> R^(d*d) such that B = G(m)

G(m) will be represented as a polynomial containing terms of the form:
I, mmT, |m| mmT, |m|^2 mmT, |m|^3 mmT ... and
|m|I, |m|^2 I, |m|^3 I, ...

where mmT denotes (m * m.T). Notice that (mmT)^2 = |m|^2 mmT, so the above terms capture everything.

Author: Arshed Nabeel
Copyright: TEE-Lab, IISc
License: MIT License

"""

import warnings

import pydaddy
from sklearn.linear_model import ridge_regression

from tensorpoly import *


class TensorFitBase(object):

    def __init__(self, threshold, order, alpha=0, smooth=True):
        """
        Initialize a fitter object

        Parameters:
             threshold: Sparsity threshold
             order: Maximum degree of the polynomial to be fit
             alpha: Regularization parameter for ridge regression
             smooth: If smooth=True, non-smooth terms (odd powers of |m|, which are not differentiable
                at m=0) will be excluded from the fit.
        """

        self.degree = order
        self.threshold = threshold
        self.alpha = alpha
        self.smooth = smooth

    def fit(self, m, y):
        """
        Fit a (vector or matrix) polynomial p over m such that p(m) = y.

        Parameters:
            m: Independent variable, a (N, d) array where d is the dimension and N is the number of data points.
            y: A (N, d) (for vector target) or (N, d, d) (for matrix target) array.

        """

        # All the corresponding arrays should be flattened out before regression.

        coeffs = self._get_coeffs()
        dictionary = self._get_dictionary(m).reshape(-1, coeffs.size)
        if self.smooth:
            dictionary = dictionary[:, self._get_smooth_terms()]
            coeffs = coeffs[self._get_smooth_terms()]

        y = y.T.ravel()

        keep = np.ones_like(coeffs, dtype=bool)
        maxiter = dictionary.shape[1]
        for it in range(maxiter):
            if np.sum(keep) == 0:
                warnings.warn('Sparsity threshold is too big, eliminated all parameters.')
                break

            coeffs_ = ridge_regression(dictionary[:, keep], y, alpha=self.alpha)
            coeffs[keep] = coeffs_
            # print(f'It: {it}, fit: {self._get_poly_object(coeffs)}')
            keep = (np.abs(coeffs) > self.threshold)
            coeffs[~keep] = 0

        if self.smooth:
            coeffs_ = self._get_coeffs()
            coeffs_[self._get_smooth_terms()] = coeffs
            coeffs = coeffs_

        return self._get_poly_object(coeffs)

    def _get_dictionary(self, m):
        raise NotImplementedError

    def _get_coeffs(self):
        raise NotImplementedError

    def _get_poly_object(self, coeffs):
        raise NotImplementedError

    def _get_non_smooth_terms(self):
        """ Terms with odd powers of |m| are not smooth (not differentiable at m=0).
            This function finds the indices of such terms. """
        raise NotImplementedError

    def _get_smooth_terms(self):
        raise NotImplementedError


class VectorFit(TensorFitBase):
    def _get_dictionary(self, m):
        return np.array([np.ones_like(m)] + [(norm(m, axis=0) ** k) * m for k in range(self.degree)]).T

    def _get_coeffs(self):
        return np.zeros(self.degree + 1)

    def _get_poly_object(self, coeffs):
        return VectorPoly(degree=self.degree, coeffs=coeffs)

    def _get_non_smooth_terms(self):
        return range(2, self.degree + 1, 2)

    def _get_smooth_terms(self):
        return range(1, self.degree + 1, 2)


class MatrixFit(TensorFitBase):
    def _get_dictionary(self, m):
        # m = np.array(m, ndmin=2)
        modm = norm(m, axis=0)
        mmT = np.einsum('ik,jk->ijk', m, m)
        I = np.eye(m.shape[0])

        i_terms = [np.einsum('i,jk->jki', modm ** k, I)
                   for k in range(self.degree + 1)]
        m_terms = [(modm ** k) * mmT
                   for k in range(self.degree - 1)]

        # The + below is a list concatenation.
        return np.array(i_terms + m_terms).T

    def _get_coeffs(self):
        return np.zeros(2 * self.degree)

    def _get_poly_object(self, coeffs):
        return MatrixPoly(degree=self.degree, coeffs=coeffs)

    def _get_non_smooth_terms(self):
        return list(range(1, self.degree + 1, 2)) + \
               list(range(self.degree + 2, 2 * self.degree, 2))

    def _get_smooth_terms(self):
        return list(range(0, self.degree + 1, 2)) + \
               list(range(self.degree + 1, 2 * self.degree, 2))


if __name__ == '__main__':
    infile = '/Users/nabeel/Data/vicsek-sde/vic_n15_eta_pi4.csv'
    data = np.loadtxt(infile, delimiter='\t')
    mx, my = data[:, 0], data[:, 1]
    # mx, my = symmetry_augmentation(mx, my, angles=8)

    subsample = 1
    dd = pydaddy.Characterize([mx, my], t=1/subsample, Dt=subsample, dt=subsample, show_summary=False)
    fx, fy = dd._ddsde._driftX_, dd._ddsde._driftY_
    gxx, gyy, gxy = dd._ddsde._diffusionX_, dd._ddsde._diffusionY_, dd._ddsde._diffusionXY_
    m = np.array([mx[:-subsample], my[:-subsample]])  # Note: This should be -Dt.
    F = np.array([fx, fy])
    G = np.array([[gxx, gxy], [gxy, gxx]])
    fitter = MatrixFit(max_degree=4, threshold=0.001, smooth=True)
    f = fitter.fit(m, G)
    print(f)
    # print(np.array(f))
