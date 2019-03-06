import os
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt




def collectCoeffs(pressure, temperature):

    # pressure in GPa
    # temperature in degK

    coeffs = {
        'alpha': 0,
          'beta': 0,
          'chi': 0,
          'delta': 0,
          'epsilon': 0
    }

    if 0.0 <= pressure <= 2:
        coeffs['alpha'] = 1.11
        coeffs['beta'] = -1.18
        coeffs['chi'] = -0.85
        coeffs['delta'] = 1680
        coeffs['epsilon'] = 487
    elif pressure == 2:
        if 2300 < temperature < 2600:
            coeffs['alpha'] = 0.84
            coeffs['beta'] = -1.22
            coeffs['chi'] = -0.85
            coeffs['delta'] = 3245
            coeffs['epsilon'] = 487
    elif pressure == 6:
        if 2300 < temperature < 2700:
            coeffs['alpha'] = 1.17
            coeffs['beta'] = -1.06
            coeffs['chi'] = -0.90
            coeffs['delta'] = 3337
            coeffs['epsilon'] = 487
    elif 2 < pressure:
        coeffs['alpha'] = 1.05
        coeffs['beta'] = -1.10
        coeffs['chi'] = -0.84
        coeffs['delta'] = 3588
        coeffs['epsilon'] = -102

    if 0.5 <= pressure <= 2:
        if 2100 < temperature < 2600:
            coeffs['alpha'] = 1.11
            coeffs['beta'] = -1.18
            coeffs['chi'] = -0.85
            coeffs['delta'] = 1680
            coeffs['epsilon'] = 487
    elif pressure == 2:
        if 2300 < temperature < 2600:
            coeffs['alpha'] = 0.84
            coeffs['beta'] = -1.22
            coeffs['chi'] = -0.85
            coeffs['delta'] = 3245
            coeffs['epsilon'] = 487
    elif pressure == 6:
        if 2300 < temperature < 2700:
            coeffs['alpha'] = 1.17
            coeffs['beta'] = -1.06
            coeffs['chi'] = -0.90
            coeffs['delta'] = 3337
            coeffs['epsilon'] = 487
    elif 2 < pressure < 18:
        if 2300 < temperature < 2700:
            coeffs['alpha'] = 1.05
            coeffs['beta'] = -1.10
            coeffs['chi'] = -0.84
            coeffs['delta'] = 3588
            coeffs['epsilon'] = -102

    return coeffs



def partition(pressure, temperature, deltaIW):
    nbo_t = 2.6
    coeffs = collectCoeffs(pressure=pressure, temperature=temperature)
    alpha = coeffs['alpha']
    beta = coeffs['beta']
    chi = coeffs['chi']
    delta = coeffs['delta']
    epsilon = coeffs['epsilon']
    logD = alpha + (beta * deltaIW) + (chi * nbo_t) + (delta * (1/temperature)) + (epsilon * (pressure/temperature))
    return logD


