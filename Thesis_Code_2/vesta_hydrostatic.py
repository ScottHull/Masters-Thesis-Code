import os
from math import sqrt, pi, exp
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt


def hydrostaticPressure(densityMelt, surfaceGrav, depth, atmosPressure):
    p = (densityMelt * surfaceGrav * depth * (10**(-9))) + atmosPressure # GPa
    return p

def adiabatTemp(Tsurf, surfaceGravity, thermExpansion, heatCapacity):
    dTdh = (thermExpansion * surfaceGravity * Tsurf) / heatCapacity
    return dTdh


magmaOceanDepth = (250) * 1000
magmaOceanDepthRange = list(np.arange(1, magmaOceanDepth, 100))
magmaOceanDepthRangeKM = [i / 1000 for i in magmaOceanDepthRange]

surfaceGrav = 0.2765738752 # m/s^2
densityMelt = 3750 # kg m3
atmosPressure = 0# GPa
Tsurf = 2000 # degK
thermExpansion = (6 * 10**(-5)) # K^-1
cp = 10**3 # J kg^-1 K^-1

pressureList = []
adiabatList = []

depthT = 1000
for index, inc in enumerate(magmaOceanDepthRange):
    if index != 0:
        dTdh = adiabatTemp(Tsurf=Tsurf, surfaceGravity=surfaceGrav, thermExpansion=thermExpansion, heatCapacity=cp) * (inc - magmaOceanDepthRange[index - 1])
        depthT += dTdh
    pressure = hydrostaticPressure(densityMelt=densityMelt, surfaceGrav=surfaceGrav, depth=inc, atmosPressure=atmosPressure)

    adiabatList.append(depthT)
    pressureList.append(pressure)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(magmaOceanDepthRangeKM, pressureList)
ax.grid()
ax.set_xlabel("Depth (km)")
ax.set_ylabel("Pressure (GPa)")
ax.set_title("Hydrostatic Pressure with Depth Below a Magma Ocean")

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(adiabatList, magmaOceanDepthRangeKM)
ax2.grid()
ax2.invert_yaxis()
ax2.set_ylabel("Depth (km)")
ax2.set_xlabel("Temperature (degK)")
ax2.set_title("Adiabatic Temperature Profile")

plt.show()