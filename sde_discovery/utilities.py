""" Visualization utilities to visualize drift and diffuion fields. """

import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
from matplotlib.ticker import ScalarFormatter
import numpy as np
from scipy.linalg import sqrtm
from scipy.signal import correlate

from tensorfit import *

CMAP = 'inferno'


def fit_model(m, drift_params, diff_params, dt=1, subsample=1, bin_avg=False):
    """
    Fit an SDE

    Args:
        m: Array containing time series data (size (2, T))
        drift_params: Dict containing parameters for drift fitting (order, threshold)
        diff_params: Dict containing parameters for diffusion fitting (order, threshold)

    Returns:
        f, g: Fitted drift and diffusion functions.
    """

    F = (m[:, subsample:] - m[:, :-subsample]) / (dt * subsample) # Instantaneous drift coefficients
    f = VectorFit(**drift_params).fit(m[:, :-subsample], F)

    r = m[:, subsample:] - m[:, :-subsample] - f(m[:, :-subsample]) * dt * subsample  # Residuals
    G = np.einsum('i...,j...->ij...', r, r) / (dt * subsample) # Instantaneous diffusion coefficients
    g = MatrixFit(**diff_params).fit(m[:, :-subsample], G)

    if bin_avg:
        # Compute binwise estimates.
        mx, my = np.linspace(-1, 1, 21), np.linspace(-1, 1, 21)
        dx, dy = mx[1] - mx[0], my[1] - my[0]

        f_bin = np.empty((2, mx.size, my.size))
        g_bin = np.empty((2, 2, mx.size, my.size))
        f_bin[:], g_bin[:] = np.nan, np.nan

        for i in range(mx.size - 1):
            for j in range(my.size - 1):
                bin_mask = ((mx[i] - dx / 2 <= m[0, :-subsample]) & (m[0, :-subsample] < mx[i] + dx / 2) &
                            (my[j] - dx / 2 <= m[1, :-subsample]) & (m[1, :-subsample] < my[j] + dx / 2))
                f_bin[:, j, i] = F[:, bin_mask].mean(axis=-1)
                g_bin[:, :, j, i] = G[:, :, bin_mask].mean(axis=-1)

        return f, g, mx, my, f_bin, g_bin

    return f, g


def plot_drift_field(ax, fx, fy, polar=True):
    """
    Visualize the drift field as a vector field.

    Args:
        ax: Axis on to which to plot to.
        fx: The (callable) drift function, that map mx, my -> f(mx, my).
        polar: Show only with i
    """

    #     r, theta = np.meshgrid(np.linspace(0, 1, 1), np.linspace(0, 2 * np.pi, 90))
    #     x, y = r * np.cos(theta), r * np.sin(theta)
    x_, y_ = np.meshgrid(np.linspace(-1, 1, 15), np.linspace(-1, 1, 15))
    xq = x_[x_ ** 2 + y_ ** 2 <= 1]
    yq = y_[x_ ** 2 + y_ ** 2 <= 1]

    xc, yc = np.meshgrid(np.linspace(-1, 1, 100), np.linspace(-1, 1, 100))
    modf = np.sqrt(fx(xc, yc) ** 2 + fy(xc, yc) ** 2)
    levels = np.linspace(0, modf[xc ** 2 + yc ** 2 <= 1].max(), 100)
    cticks = levels[::33]
    cf = ax.contourf(xc, yc, modf, levels=levels, cmap=CMAP)
    for c in cf.collections:
        c.set_rasterized(True)

    qv = ax.quiver(xq, yq, fx(xq, yq), fy(xq, yq), color=[1, 1, 1, 1],
                   width=0.008, pivot='tail')
    ax.set(xlabel='$m_x$', ylabel='$m_y$',)# title='Drift Field')
    # ax.set(xticks=[-1, -.5, 0, .5, 1], yticks=[-1, -.5, 0, .5, 1])
    ax.set_aspect('equal', 'box')

    cax = plt.colorbar(cf, ax=ax, fraction=0.0453, ticks=cticks,
                       format=ScalarFormatter(useMathText=True))
    cax.ax.ticklabel_format(style='scientific', scilimits=(0, 0), useMathText=True)


def plot_diffusion_field(ax, gxx, gyy, gxy, scale=1):
    """
    Visualize the diffusion field

    Args:
        ax: Axis on to which to plot to
        gxx, gyy: Diffusion functions
        gxy: Cross diffusion
        scale: Scale factor to scale the ellipses with
    """

    x_, y_ = np.meshgrid(np.linspace(-1, 1, 15), np.linspace(-1, 1, 15))
    xs = x_[x_ ** 2 + y_ ** 2 <= 1]
    ys = y_[x_ ** 2 + y_ ** 2 <= 1]
    xy = np.column_stack((xs.ravel(), ys.ravel()))
    maj_axis = np.empty(xs.size)
    min_axis = np.empty(xs.size)
    angles = np.empty(xs.size)
    colors = np.empty(xs.size)
    for i, (x, y) in enumerate(zip(xs.ravel(), ys.ravel())):
        diff = [[gxx(x, y), gxy(x, y)],
                [gxy(x, y), gyy(x, y)]]

        # Eigendecomposition is computed using the lower triangular part, assuming symmetry
        eigval, eigvec = np.linalg.eigh(diff, UPLO='L')
        maj_axis[i] = scale * eigval[0]
        min_axis[i] = scale * eigval[1]
        angles[i] = np.arctan2(eigvec[0, 1], eigvec[0, 0]) * 180 / np.pi
        colors[i] = np.prod(eigval)

    ec = EllipseCollection(maj_axis, min_axis, angles,
                           offsets=xy, offset_transform=ax.transData,
                           linewidths=1.5,
                           # edgecolors=[1, 1, 1, 1],
                           # edgecolors=plt.cm.inferno(plt.Normalize()(colors)),
                           facecolors=(0, 0, 0, 0),
                           cmap=CMAP)
    ec.set_array(colors)
    ax.add_collection(ec)
    ax.autoscale_view()
    # ax.scatter(xs.ravel(), ys.ravel(), marker='+', color=(0.8, 0.8, 0.8))
    ax.set(xlabel='$m_x$', ylabel='$m_y$', )# title='Diffusion Field')
    ax.set(xticks=[-1, -.5, 0, .5, 1], yticks=[-1, -.5, 0, .5, 1])
    ax.set_aspect('equal', 'box')

    cax = plt.colorbar(ec, ax=ax, fraction=0.0453)
    cax.ax.ticklabel_format(style='scientific', scilimits=(-3, 3), useMathText=True)


def symmetry_augmentation(mx, my, angles=4, mirror=False):

    if mirror:
        mx = np.concatenate(mx, -mx, mx)
        my = np.concatenate(my, my, -my)

    mx_, my_ = np.empty(mx.size * angles), np.empty(my.size * angles)
    mx_[:mx.size] = mx
    my_[:my.size] = my

    phis = np.linspace(2 * np.pi / angles, (1 - 1 / angles) * 2 * np.pi, angles - 1)
    for i, phi in enumerate(phis):
        r = np.sqrt(mx ** 2 + my ** 2)
        theta = np.arctan2(my, mx)

        mx_[(mx.size * (i + 1)):(mx.size * (i + 2))] = r * np.cos(theta - phi)
        my_[(my.size * (i + 1)):(my.size * (i + 2))] = r * np.sin(theta - phi)

    return mx_, my_


def simulate(f, G, dt, length, x0=0.):
    """ Generate a simulated time series using an SDE with drift f and diffusion sqrt(G). """

    t = np.linspace(0, length*dt, length)
    x = np.zeros((2, length))
    x[:, 0] = x0
    # dW = np.random.normal(size=(2, length))

    for i in range(length - 1):
        df = dt * f(x[:, i])
        dg = np.sqrt(dt) * sqrtm(G(x[:, i])) @ np.random.normal(size=(2, ))
        x_next = x[:, i] + df + dg
        if np.linalg.norm(x_next) <= 1:
            x[:, i + 1] = x_next
        else:
            x[:, i + 1] = x_next / np.linalg.norm(x_next)

        # x[:, i + 1] = x_next #if (x_next[0] ** 2 + x_next[1] ** 2 <= 1) else 1

    return t, x


def compute_acf(x):
    """
    Calculates autocorrelation using wiener khinchin theorem.
    """

    x = x - x.mean()
    c = correlate(x, x)
    c = c[c.shape[0] // 2:]
    c /= c[0]
    return c


def compute_r2(f, x, y):
    """
    Compute the R^2 value for a fitted function f, based on the given values fx
    evaluated at points x.
    
    Args:
        f : Fitted function
        x (np.ndarray): Array of x-values to evaluate the function at.
        y (np.ndarray): Array of same size as x, with 'actual' values of fx
        
    Returns:
        (float) R^2 value.
    """
    
    r2 = 1 - np.nanmean((y - f(x)) ** 2, axis=-1) / np.nanvar(y, axis=-1)
    return np.nanmean(r2)
    
    
    
    