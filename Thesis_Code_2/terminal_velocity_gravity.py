import numpy as np
from math import sqrt, pi
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt

def velocity(gravAccel, diameter, densityDroplet, densityMelt):
    Cd = 0.2
    v = sqrt(((4 * gravAccel * diameter) / (3 * Cd)) * ((densityDroplet - densityMelt) / densityMelt))
    return v


gravity_earth = 9.8
gravity_vesta = 0.25

gravity_range = np.arange(0.1, 10.1, 0.1)
density_melt_range = np.arange(2000, 5000, 100)



dropletRadius = (1 / 100) / 2 # 0.5cm
dropletDiameter = dropletRadius * 2
volumeDroplet = (4/3) * pi * (dropletRadius**3)
densityDroplet = 7800 # kg m3
massDroplet = densityDroplet * volumeDroplet

viscosityMelt = 10**(-2)
densityMelt = 3750  # kg m3

results_gravity = []
results_earth_density = []
results_vesta_density = []

if __name__ == "__main__":

    for i in gravity_range:
        v = velocity(gravAccel=i, diameter=dropletDiameter, densityMelt=densityMelt, densityDroplet=densityDroplet)
        results_gravity.append(v)
    for i in density_melt_range:
        v_e = velocity(gravAccel=gravity_earth, diameter=dropletDiameter, densityMelt=i, densityDroplet=densityDroplet)
        v_v = velocity(gravAccel=gravity_vesta, diameter=dropletDiameter, densityMelt=i, densityDroplet=densityDroplet)
        results_earth_density.append(v_e)
        results_vesta_density.append(v_v)

    fig_density = plt.figure()
    fig_gravity = plt.figure()
    ax_density = fig_density.add_subplot(111)
    ax_gravity = fig_gravity.add_subplot(111)

    ax_gravity.plot(gravity_range, results_gravity)

    ax_density.plot(density_melt_range, results_earth_density, label="g = 9.8 m/s")
    ax_density.plot(density_melt_range, results_vesta_density, label="g = 0.25 m/s")

    ax_gravity.grid()
    ax_density.grid()

    ax_density.legend()

    plt.show()