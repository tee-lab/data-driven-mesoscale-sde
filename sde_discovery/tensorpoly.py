"""
Utilities for expressing vector and matrix polynomial functions generated from a vector m.

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

import numpy as np
from numpy.linalg import norm


class TensorPolyBase(object):
    def __init__(self, degree, coeffs):
        self.coeffs = np.array(coeffs)
        self.degree = degree

    def __str__(self):
        return self.__repr__()

    def __array__(self):
        return self.coeffs

    def __len__(self):
        return self.coeffs.size

    def __repr__(self):
        raise NotImplementedError

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class VectorPoly(TensorPolyBase):

    def __init__(self, degree, coeffs):
        assert len(coeffs) == degree + 1, \
            f'For degree {degree}, number of coefficients should be {degree + 1}.'

        super().__init__(degree, coeffs)

    def __repr__(self):
        def term(n):
            if n == 0:
                return ''
            elif n == 1:
                return 'm'
            elif n == 2:
                return f'|m| m'
            else:
                return f'|m|^{n - 1} m'

        terms = [term(n) for n in range(self.degree + 1)]
        terms_with_coeffs = [f'{c:.3f}{t}' for (c, t) in zip(self.coeffs, terms) if c != 0]

        if terms_with_coeffs:
            return ' + '.join(terms_with_coeffs)
        else:
            return '0'

    def __call__(self, m):
        m = np.array(m)
        modm = norm(m, axis=0)
        terms = np.array([np.ones_like(m)] + [(modm ** k) * m for k in range(self.degree)])
        terms_with_coeffs = terms.T * self.coeffs
        return terms_with_coeffs.sum(axis=-1).T


class MatrixPoly(TensorPolyBase):

    def __init__(self, degree, coeffs):
        assert len(coeffs) == 2 * degree, \
            f'For degree {degree}, number of coefficients should be {2 * degree}.'

        super().__init__(degree, coeffs)

    def __repr__(self):

        def i_term(n):
            if n == 0:
                return 'I'
            elif n == 1:
                return '|m| I'
            else:
                return f'|m|^{n} I'

        def m_term(n):
            assert n >= 2

            if n == 2:
                return 'mmᵀ'
            elif n == 3:
                return '|m| mmᵀ'
            else:
                return f'|m|^{n - 2} mmᵀ'

        terms = [i_term(k) for k in range(self.degree + 1)] + \
                [m_term(k) for k in range(2, self.degree + 1)]

        terms_with_coeffs = [f'{c:.3f}{t}' for (c, t) in zip(self.coeffs, terms) if c != 0]

        if terms_with_coeffs:
            return ' + '.join(terms_with_coeffs)
        else:
            return '0'

    def __call__(self, m):
        # m = np.expand_dims(np.array(m), 0)
        m = np.array(m)
        modm = norm(m, axis=0)
        mmT = np.einsum('i...,j...->ij...', m, m)
        I = np.eye(m.shape[0])

        i_terms = [np.einsum('...,jk->jk...', modm ** k, I) for k in range(self.degree + 1)]
        m_terms = [(modm ** k) * mmT for k in range(self.degree - 1)]

        # The + below is a list concatenation.
        terms_with_coeffs = self.coeffs * np.array(i_terms + m_terms).T
        return terms_with_coeffs.sum(axis=-1).T


if __name__ == '__main__':
    # poly = VectorPoly(degree=5, coeffs=[1, 0, 0, 0, 0, 0])
    # print(poly([1, 2, 4]))
    # print(poly)
    G = MatrixPoly(degree=3, coeffs=[0, 0, 0, 0, 1, 0, ])
    print(G)
    x, y = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 10))
    print(G([1, 2]))
    print(G([x, y]))
