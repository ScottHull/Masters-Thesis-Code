import os
from math import sqrt, pi, exp
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt


def accelerationGravity(radius, meltDensity):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    earthMass = 5.972 * 10**(24) # kg
    # Rearth = 6371 * 1000 # m
    Rearth = (250) * 1000 # m
    gSurface = (4/3) * pi * G * meltDensity * Rearth
    gAccel = gSurface * (1 - (radius/Rearth))
    return gAccel

def forceBuoyant(gravAccel, diameter, meltDensity):
    f = (1/6) * (meltDensity * pi * (diameter**3) * gravAccel)
    return f

def forceGravity(gravAccel, dropletDensity, diameter):
    f = (1/6) * (dropletDensity * pi * (diameter**3) * gravAccel)
    return f

def forceDrag(diameter, densityMelt, terminalVelocity):
    Cd = 0.2
    r = diameter / 2
    A = pi * (r**2)
    f = (1/2) * (Cd * A * densityMelt * (terminalVelocity**2))
    return f

def velocity(gravAccel, diameter, densityDroplet, densityMelt):
    Cd = 0.2
    v = sqrt(((4 * gravAccel * diameter) / (3 * Cd)) * ((densityDroplet - densityMelt) / densityMelt))
    return v

def frictionCoeff(gravAccel, densityMelt, densityDroplet, diameter, viscosityMelt):
    f = (pi/6) * ((densityDroplet - densityMelt) / densityMelt) * ((densityMelt / viscosityMelt)**2) * gravAccel * (diameter**3)
    return f

def rubiePartitioning(diameter, initialConcMetal, initialConcSilicate, velocity, densityDroplet, densityMelt, distributionCoeff, time):
    radius = diameter / 2
    diffusivityMelt = 0
    Fm = ((radius**3) * (densityDroplet / 3)) / (((radius**3) * (densityDroplet / 3)) + ((radius**2) * densityMelt) * sqrt((2 * diffusivityMelt * radius)/ velocity))
    Fs = 1 - Fm
    A = (velocity / (2 * radius)) * (Fm / (Fm + (Fs / distributionCoeff)))
    B = (velocity / (2 * radius)) * ((Fs * initialConcSilicate) / (Fm / (Fm + (Fs / distributionCoeff))))
    Cm = ((initialConcMetal + (B / A)) * exp(A * time)) - (B / A)
    return Cm

def myPartitioning(distributionCoeff, Dpred, Dcurr, conc):
    pass


dropletRadius = (1 / 100) / 2 # 0.5cm
dropletDiameter = dropletRadius * 2
volumeDroplet = (4/3) * pi * (dropletRadius**3)
densityDroplet = 7800 # kg m3
massDroplet = densityDroplet * volumeDroplet
densityMelt = 3750 # kg m3
viscosityMelt = 10**(-2)

magmaOceanDepth = (250) * 1000
magmaOceanDepthRange = list(reversed(np.arange(1, magmaOceanDepth, 1000)))
magmaOceanDepthRangeKM = [i / 1000 for i in magmaOceanDepthRange]

fig = plt.figure()
ax = fig.add_subplot(111)

vList = []
fList = []
fbList = []
fgList = []
fdList = []




for inc in magmaOceanDepthRange:
    g = accelerationGravity(radius=inc, meltDensity=densityMelt)
    f = frictionCoeff(gravAccel=g, densityMelt=densityMelt, densityDroplet=densityDroplet, diameter=dropletDiameter, viscosityMelt=viscosityMelt)
    v = velocity(gravAccel=g, diameter=dropletDiameter, densityDroplet=densityDroplet, densityMelt=densityMelt)
    fb = forceBuoyant(gravAccel=g, diameter=dropletDiameter, meltDensity=densityMelt)
    fg = forceGravity(gravAccel=g, dropletDensity=densityDroplet, diameter=dropletDiameter)
    fd = forceDrag(diameter=dropletDiameter, densityMelt=densityMelt, terminalVelocity=v)
    balance = (fd + fb) - fg

    vList.append(v)
    fList.append(f)
    fbList.append(fb)
    fgList.append(fg)
    fdList.append(fd)


ax.plot(magmaOceanDepthRangeKM, fbList, label="Buoyant Force", c='b', linestyle='--')
ax.plot(magmaOceanDepthRangeKM, fgList, label="Gravitational Force", c='g', linestyle='--')
ax.plot(magmaOceanDepthRangeKM, fdList, label="Drag Force", c='r', linestyle='--')

ax2 = ax.twinx()
ax2.plot(magmaOceanDepthRangeKM, vList, label="Droplet Terminal Velocity", c='black')
ax2.set_ylabel("Terminal Velocity (m/s)")

ax.grid()
fig.legend(loc='upper right', bbox_to_anchor=(0.9,0.85))
ax.set_xlabel("Depth Below Magma Ocean Surface (km)")
ax.set_ylabel("Force (N)")
ax.set_title("Forces Acting On An Iron Droplet At Terminal Velocity")
plt.rcParams.update({'font.size': 80})




plt.show()


