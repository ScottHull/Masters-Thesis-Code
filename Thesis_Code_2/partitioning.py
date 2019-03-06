import os
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt
import pandas as pd

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
    # elif pressure == 2:
    #     if 2300 < temperature < 2600:
    #         coeffs['alpha'] = 0.84
    #         coeffs['beta'] = -1.22
    #         coeffs['chi'] = -0.85
    #         coeffs['delta'] = 3245
    #         coeffs['epsilon'] = 487
    # elif pressure == 6:
    #     if 2300 < temperature < 2700:
    #         coeffs['alpha'] = 1.17
    #         coeffs['beta'] = -1.06
    #         coeffs['chi'] = -0.90
    #         coeffs['delta'] = 3337
    #         coeffs['epsilon'] = 487
    elif 2 < pressure:
        coeffs['alpha'] = 1.05
        coeffs['beta'] = -1.10
        coeffs['chi'] = -0.84
        coeffs['delta'] = 3588
        coeffs['epsilon'] = -102

    # if 0.5 <= pressure <= 2:
    #     if 2100 < temperature < 2600:
    #         coeffs['alpha'] = 1.11
    #         coeffs['beta'] = -1.18
    #         coeffs['chi'] = -0.85
    #         coeffs['delta'] = 1680
    #         coeffs['epsilon'] = 487
    # # elif pressure == 2:
    # #     if 2300 < temperature < 2600:
    # #         coeffs['alpha'] = 0.84
    # #         coeffs['beta'] = -1.22
    # #         coeffs['chi'] = -0.85
    # #         coeffs['delta'] = 3245
    # #         coeffs['epsilon'] = 487
    # # elif pressure == 6:
    # #     if 2300 < temperature < 2700:
    # #         coeffs['alpha'] = 1.17
    # #         coeffs['beta'] = -1.06
    # #         coeffs['chi'] = -0.90
    # #         coeffs['delta'] = 3337
    # #         coeffs['epsilon'] = 487
    # elif 2 < pressure < 18:
    #     if 2300 < temperature < 2700:
    #         coeffs['alpha'] = 1.05
    #         coeffs['beta'] = -1.10
    #         coeffs['chi'] = -0.84
    #         coeffs['delta'] = 3588
    #         coeffs['epsilon'] = -102

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
Tsurf = 2000 # degK
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
    logDlist = []
    for index, j in enumerate(magmaOceanDepthRange):
        if index != 0:
            logD = partition(pressure=pressureList[index], temperature=adiabatList[index], deltaIW=i)
            logDlist.append(logD)
    logfO2list.append(logDlist)

fig = plt.figure()
ax = fig.gca()
ax2 = ax.twiny()
for index, i in enumerate(logfO2list):
    ax.plot(pressureList[1:], i, label="IW {}".format(deltaIW[index]))
    ax.plot(adiabatList[1:], i, label="IW {}".format(deltaIW[index]))


    ax2.plot(magmaOceanDepthRangeKM[1:], i)
    ax2.set_ylabel('Depth')

ax.grid()
ax.legend(loc='upper right')
# ax.set_xlabel("Pressure (GPa)")
ax.set_xlabel("Temperature (degK)")
ax.set_ylabel("log(D)")
# ax2.set_xlabel("Depth (km)")

df = pd.read_csv('w.csv')

for row in df.index:
    pressure = df['Pressure'][row]
    temperature = df['Temperature'][row]
    fO2 = df['fO2']
    measuredlogD = df['Corrected_logD'][row]
    capsule = df['Capsule'][row]
    if capsule == 'C':
        # ax.scatter(pressure, measuredlogD, color='blue')
        ax.scatter(temperature, measuredlogD, color='blue')
    else:
        # ax.scatter(pressure, measuredlogD, color='red')
        ax.scatter(temperature, measuredlogD, color='red')

# ax.set_title("Partitioning Model (Cottrell et al. 2009)")

plt.show()
