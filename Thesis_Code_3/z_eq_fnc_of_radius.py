from math import pi, log, sqrt
import pandas as pd
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})



def rFromWeberTurbulent(density_melt, density_droplet, gravity, Cd=0.2, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = sqrt(weber_number * (((16 * gravity * (density_droplet - density_melt)) / (3 * Cd * surface_tension)) *
                        ((density_droplet - density_melt) / density_melt))**(-1))
    return r

def z_eq(radius_droplet, dynamic_viscosity, settling_velocity, diffusion_coeff, density_melt):
    kinematic_viscosity = dynamic_viscosity / density_melt
    num = ((2 * radius_droplet)**(3/2)) * (kinematic_viscosity**(1/6)) * (settling_velocity**(1/2))
    den = diffusion_coeff**(2/3)
    z = 4.624 * (num/den)
    return z

def turbulentVelocity(gravity, droplet_radius, density_droplet, density_melt, Cd=0.2):
    diameter = 2 * droplet_radius
    v = sqrt(((4 * gravity * diameter) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v


radii = np.arange(0.001, 0.02 + 0.001, 0.001)
gravity = 9.8
density_droplet = 7800
density_melt = 3750
eta = 10**-2
diffusion_coeff = 10**-9

velocities = []
z_eqs = []

for i in radii:
    v = turbulentVelocity(gravity=gravity, density_melt=density_melt, density_droplet=density_droplet, droplet_radius=i)
    z = z_eq(radius_droplet=i, dynamic_viscosity=eta, settling_velocity=v, diffusion_coeff=diffusion_coeff,
             density_melt=density_melt)
    velocities.append(v)
    z_eqs.append(z)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

