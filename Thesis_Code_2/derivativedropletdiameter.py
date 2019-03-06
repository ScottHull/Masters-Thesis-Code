import os
from math import sqrt, pi, exp
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt


def derivHard(gravAccel, diameter, densityDroplet, densityMelt):
    Cd = 0.2
    a = (1/2) * (((4*gravAccel*diameter)/(3*Cd)) * ((densityDroplet - densityMelt)/densityMelt))**(-1/2)
    b = ((4*gravAccel)/(3*Cd)) * ((densityDroplet - densityMelt)/densityMelt)
    dVdD = a * b
    return dVdD

def derivEasy(gravAccel, diameter, densityDroplet, densityMelt):
    Cd = 0.2
    dVdD = (sqrt(gravAccel) * sqrt(Cd * densityMelt) * sqrt((densityDroplet - densityMelt))) / (sqrt(3) * Cd * densityMelt * sqrt(diameter))
    return dVdD

def velocity(gravAccel, diameter, densityDroplet, densityMelt):
    Cd = 0.2
    v = sqrt(((4 * gravAccel * diameter) / (3 * Cd)) * ((densityDroplet - densityMelt) / densityMelt))
    return v


dropletRadius = [((i / 100) / 2) for i in np.arange(0.5, 100, 1)] # m
dropletDiameter = [i * 2 for i in dropletRadius]
volumeDroplet = [(4/3) * pi * (i**3) for i in dropletRadius]
densityDroplet = 7800 # kg m3
massDroplet = [densityDroplet * i for i in volumeDroplet]
densityMelt = 3750 # kg m3
viscosityMelt = 10**(-2)
gravAccel = 9.8

vList = []
dVdDList = []


for i in dropletDiameter:
    v = velocity(gravAccel=gravAccel, diameter=i, densityMelt=densityMelt, densityDroplet=densityDroplet)
    de = derivEasy(gravAccel=gravAccel, densityDroplet=densityDroplet, densityMelt=densityMelt, diameter=i)
    dh = derivHard(gravAccel=gravAccel, densityDroplet=densityDroplet, densityMelt=densityMelt, diameter=i)
    vList.append(v)
    dVdDList.append(de)


testslope = derivEasy(gravAccel=gravAccel, diameter=0.1, densityMelt=densityMelt, densityDroplet=densityDroplet)
testvelocity = velocity(gravAccel=gravAccel, diameter=0.1, densityMelt=densityMelt, densityDroplet=densityDroplet)
b = testvelocity - (testslope * 0.1)
tanline = [(testslope * i) + b for i in dropletDiameter]



fig = plt.figure()
ax = fig.gca()

ax.plot(dropletDiameter, vList, label='Terminal Velocity')
ax.plot(dropletDiameter, tanline, color='black', linestyle='--', linewidth=2.0)
ax2 = ax.twinx()
ax2.plot(dropletDiameter, dVdDList, linestyle='--', label='dv/dd')
ax.set_xlabel("Droplet Diameter (m)")
ax.set_ylabel("Terminal Velocity")
ax2.set_ylabel("d(velocity)/d(droplet diameter)")
ax.grid()
ax.legend()

plt.show()

