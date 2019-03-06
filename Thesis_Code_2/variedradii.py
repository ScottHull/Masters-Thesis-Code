import os
from math import sqrt, pi, exp
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt


def accelerationGravity(depth, meltDensity):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    earthMass = 5.972 * 10**(24) # kg
    # Rearth = 6371 * 1000 # m
    Rearth = (3500 + 2900) * 1000 # m
    gSurface = (4/3) * pi * G * meltDensity * Rearth
    gAccel = gSurface * (1 - (depth/Rearth)) # depth below Earth's surface
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

def derivEasy(gravAccel, diameter, densityDroplet, densityMelt):
    # change in terminal velocity with respect to droplet diameter
    Cd = 0.2
    dVdD = (sqrt(gravAccel) * sqrt(Cd * densityMelt) * sqrt((densityDroplet - densityMelt))) / (sqrt(3) * Cd * densityMelt * sqrt(diameter))
    return dVdD


dropletRadius = [((i / 100) / 2) for i in np.arange(0.5, 2.5, 0.5)] # m
dropletDiameter = [i * 2 for i in dropletRadius]
volumeDroplet = [(4/3) * pi * (i**3) for i in dropletRadius]
densityDroplet = 7800 # kg m3
massDroplet = [densityDroplet * i for i in volumeDroplet]
densityMelt = 3750 # kg m3
viscosityMelt = 10**(-2)

magmaOceanDepth = (3500 + 2900) * 1000
magmaOceanDepthRange = list(reversed(np.arange(1, magmaOceanDepth, 10000)))
magmaOceanDepthRangeKM = [i / 1000 for i in magmaOceanDepthRange]

fig = plt.figure()
ax = fig.gca()


vList = []
fList = []
fbList = []
fgList = []
fdList = []
gList = []

calcGrav = False

for index, i in enumerate(dropletRadius):

    d = dropletDiameter[index]
    v = volumeDroplet[index]
    m = massDroplet[index]

    vListr = []
    fListr = []
    fbListr = []
    fgListr = []
    fdListr = []


    for inc in magmaOceanDepthRange:
        g = accelerationGravity(depth=inc, meltDensity=densityMelt)
        f = frictionCoeff(gravAccel=g, densityMelt=densityMelt, densityDroplet=densityDroplet, diameter=d, viscosityMelt=viscosityMelt)
        v = velocity(gravAccel=g, diameter=d, densityDroplet=densityDroplet, densityMelt=densityMelt)
        fb = forceBuoyant(gravAccel=g, diameter=d, meltDensity=densityMelt)
        fg = forceGravity(gravAccel=g, dropletDensity=densityDroplet, diameter=d)
        fd = forceDrag(diameter=d, densityMelt=densityMelt, terminalVelocity=v)
        balance = (fd + fb) - fg

        if calcGrav is False:
            gList.append(g)
        vListr.append(v)
        fListr.append(f)
        fbListr.append(fb)
        fgListr.append(fg)
        fdListr.append(fd)

    calcGrav = True
    vList.append(vListr)
    fList.append(fListr)
    fbList.append(fbListr)
    fgList.append(fbListr)
    fdList.append(fbListr)


# ax.plot(magmaOceanDepthRangeKM, fbList, label="Buoyant Force", c='b', linestyle='--')
# ax.plot(magmaOceanDepthRangeKM, fgList, label="Gravitational Force", c='g', linestyle='--')
# ax.plot(magmaOceanDepthRangeKM, fdList, label="Drag Force", c='r', linestyle='--')

for index, i in enumerate(dropletRadius):
    ax.plot(vList[index], magmaOceanDepthRangeKM, label="r = {} cm".format(dropletRadius[index] * 100))
ax2 = ax.twiny()
ax2.plot(gList, magmaOceanDepthRangeKM, linestyle='--', label="grav. accel.", color='black', linewidth=3)
ax.set_ylabel("Depth Below Earth Surface (km)")
ax2.set_xlabel("Gravitational Acceleration")

ax.plot(list(np.arange(0, 2, 0.2)), tanline)

ax.invert_yaxis()

ax.grid()
fig.legend(loc='upper right')
ax.set_xlabel("Terminal Velocity (m/s)")

plt.show()