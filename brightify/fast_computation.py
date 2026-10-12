# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

#%% Fast computation functions

def to_f(x):
    """
    Convert Pandas DataFrame or Series to a NumPy array.
    If input is already a NumPy array, return as is.
    """
    if isinstance(x, (pd.DataFrame, pd.Series)):
        return x.to_numpy()  
    return x

def multiply_f(x, y):
    """
    Element-wise multiplication.
    If y is a scalar, multiply directly.
    If y is an array, expand dimensions for broadcasting.
    """
    if y.ndim == 0:
        return to_f(x) * y
    return to_f(x) * to_f(y)[:, np.newaxis]

def divide_f(x, y):
    """
    Element-wise division.
    If y is a scalar, divide directly.
    If y is an array, expand dimensions for broadcasting.
    """
    if y.ndim == 0:
        return to_f(x) / y
    return to_f(x) / to_f(y)[:, np.newaxis]

def sum_f(x):
    """
    Compute the sum along axis 0.
    """
    return np.sum(to_f(x), axis=0)

def mean_f(x, w):
    """
    Compute the weighted mean of x using weights w.
    """
    return divide_f(sum_f(multiply_f(x, w)), sum_f(w))

def norm_f(x):
    """
    Compute the Euclidean norm (magnitude) along the last axis.
    """
    ax = x.ndim - 1
    return np.linalg.norm(to_f(x), axis=ax)

def direction_f(x):
    """
    Compute the unit vector (direction) of x.
    """
    return divide_f(x, norm_f(x))

def dot_f(x, y):
    """
    Compute the dot product of two arrays.
    """
    return np.dot(to_f(x), to_f(y))

def tangent_vector_f(x, y):
    """
    Compute the tangent vector projection of x onto y.
    """
    return np.outer(dot_f(x, y), to_f(y))

def normal_vector_f(x, y):
    """
    Compute the normal (perpendicular) component of x relative to y.
    """
    return to_f(x) - tangent_vector_f(x, y)

def cos_f(x, y):
    """
    Compute the cosine of the angle between two vectors.
    """
    return dot_f(x, y) / (norm_f(x) * norm_f(y))

def arccos_f(x):
    """
    Compute the inverse cosine (arccos) of x.
    """
    return np.arccos(to_f(x))

def spher2cart(r, theta, phi):
    """
    Convert spherical coordinates (r, theta, phi)
    to Cartesian coordinates (x, y, z).
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x, y, z]).T

def cart2spher(x, y, z):
    """
    Convert Cartesian coordinates (x, y, z) to
    spherical coordinates (r, theta, phi).
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    return  np.array([r, theta, phi]).T

def cosRR_f(x, y):
    """
    Compute the cosine similarity for multiple vector pairs.
    This function processes each row of `x` and `y` as separate vectors
    and computes the cosine similarity in a batch-wise manner.
    """
    return np.sum(to_f(x) * to_f(y), axis=1) / (norm_f(x) * norm_f(y))