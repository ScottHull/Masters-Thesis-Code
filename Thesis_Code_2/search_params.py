import numpy as np
from math import pi, sqrt, log
import matplotlib.pyplot as plt



def adiabatTemp(Tsurf, surfaceGravity, thermExpansion, heatCapacity):
    dTdh = (thermExpansion * surfaceGravity * Tsurf) / heatCapacity
    return dTdh

def hydrostaticPressure(densityMelt, surfaceGrav, depth, atmosPressure):
    p = (densityMelt * surfaceGrav * depth * (10**(-9))) + atmosPressure # GPa
    return p

def calcSurfaceGravity(body_mass, body_radius):
    grav_constant = 6.674 * (10**(-11))
    g = (grav_constant * body_mass) / (body_radius**2)
    return g

def calcDepthGravity(body_radius, depth_radius, body_surf_grav):
    gAccel = body_surf_grav * (1 - (depth_radius/body_radius))
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




dropletRadius = (1 / 100) / 2 # 0.5cm
dropletDiameter = dropletRadius * 2
volumeDroplet = (4/3) * pi * (dropletRadius**3)
densityDroplet = 7800 # kg m3
massDroplet = densityDroplet * volumeDroplet


mass_vesta = 2.589 * (10**20) # kg
radius_vesta = 265 * 1000 # m
surfaceGrav = calcSurfaceGravity(body_mass=mass_vesta, body_radius=radius_vesta)  # m/s^2

viscosityMelt = 10**(-2)
densityMelt = 3750  # kg m3
# atmosPressure = (1.01325 * (10 ** (-4)))  # GPa
Tsurf = 2000  # degK
thermExpansion = (6 * 10 ** (-5))  # K^-1
cp = 10 ** 3  # J kg^-1 K^-1
deltaIW = [0, -1, -2]

diffusion = 10**(-8)

pressureList = []
adiabatList = []
logfO2List = []
velocityList = []

depthT = Tsurf

magmaOceanDepth = radius_vesta
magmaOceanDepthRange = list(np.arange(1, magmaOceanDepth, 1000))
magmaOceanDepthRangeKM = [i / 1000 for i in magmaOceanDepthRange]



if __name__ == "__main__":

    vesta_surface_g = calcSurfaceGravity(body_radius=radius_vesta, body_mass=mass_vesta)

    for index, inc in enumerate(magmaOceanDepthRange):
        g = calcDepthGravity(body_surf_grav=vesta_surface_g, body_radius=radius_vesta, depth_radius=radius_vesta - inc)
        v = velocity(densityDroplet=densityDroplet, diameter=dropletDiameter, densityMelt=densityMelt, gravAccel=g)
        velocityList.append(sqrt(diffusion * dropletDiameter / v))


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(magmaOceanDepthRangeKM, velocityList)
    ax.grid()
    ax.set_title("Characteristic Chemical Reaction for a 1cm Iron Droplet")
    ax.set_xlabel("Depth (km)")
    ax.set_ylabel("Reaction time (s)")
    plt.gca().invert_xaxis()
    plt.show()

#
#     for index, inc in enumerate(magmaOceanDepthRange):
#         if index != 0:
#             dTdh = adiabatTemp(Tsurf=Tsurf, surfaceGravity=surfaceGrav, thermExpansion=thermExpansion,
#                                heatCapacity=cp) * (inc - magmaOceanDepthRange[index - 1])
#             depthT += dTdh
#         pressure = hydrostaticPressure(densityMelt=densityMelt, surfaceGrav=surfaceGrav, depth=inc,
#                                        atmosPressure=atmosPressure)
#
#         adiabatList.append(depthT)
#         pressureList.append(pressure)
#
#     print(adiabatList)