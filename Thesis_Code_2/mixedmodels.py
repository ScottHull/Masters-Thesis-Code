import os
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt


def cottrellModel(pressure, temperature, deltaIW, nbo_t=2.6):

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

    alpha = coeffs['alpha']
    beta = coeffs['beta']
    chi = coeffs['chi']
    delta = coeffs['delta']
    epsilon = coeffs['epsilon']
    logD = alpha + (beta * deltaIW) + (chi * nbo_t) + (delta * (1 / temperature)) + (epsilon * (pressure / temperature))

    return logD


def seibertModel(pressure, temperature, deltaIW, nbo_t=2.6):

    alpha = 1.96
    beta = -937
    chi = -55
    delta = -0.57
    epislon = 0

    logD = alpha + (beta* (1 / temperature)) + (chi * (pressure / temperature))

    return logD



def hydrostaticPressure(densityMelt, surfaceGrav, depth, atmosPressure):
    p = (densityMelt * surfaceGrav * depth * (10**(-9))) + atmosPressure # GPa
    return p

def adiabatTemp(Tsurf, surfaceGravity, thermExpansion, heatCapacity):
    dTdh = (thermExpansion * surfaceGravity * Tsurf) / heatCapacity
    return dTdh



magmaOceanDepth = (3500 + 2900) * 1000
magmaOceanDepthRange = list(np.arange(1, magmaOceanDepth, 10000))
magmaOceanDepthRangeKM = [i / 1000 for i in magmaOceanDepthRange]

surfaceGrav = 9.8 # m/s^2
densityMelt = 3750 # kg m3
atmosPressure = (1.01325 * (10**(-4))) # GPa
Tsurf = 1000 # degK
thermExpansion = (6 * 10**(-5)) # K^-1
cp = 10**3 # J kg^-1 K^-1
deltaIW = [0, -1, -2]

pressureList = []
adiabatList = []
logfO2list = []

depthT = 1000
for index, inc in enumerate(magmaOceanDepthRange):
    if index != 0:
        dTdh = adiabatTemp(Tsurf=Tsurf, surfaceGravity=surfaceGrav, thermExpansion=thermExpansion, heatCapacity=cp) * (inc - magmaOceanDepthRange[index - 1])
        depthT += dTdh
    pressure = hydrostaticPressure(densityMelt=densityMelt, surfaceGrav=surfaceGrav, depth=inc, atmosPressure=atmosPressure)

    adiabatList.append(depthT)
    pressureList.append(pressure)

for i in deltaIW:
    cottrell_logDlist = []
    seibert_logDlist = []
    for index, j in enumerate(magmaOceanDepthRange):
        if index != 0:
            cottrell_logD = cottrellModel(pressure=pressureList[index], temperature=adiabatList[index], deltaIW=i)
            seibert_logD = seibertModel(pressure=pressureList[index], temperature=adiabatList[index], deltaIW=i)
            cottrell_logDlist.append(cottrell_logD)
            seibert_logDlist.append(seibert_logD)
    logfO2list.append([cottrell_logDlist, seibert_logDlist])

