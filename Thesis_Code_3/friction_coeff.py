import os
from math import sqrt, pi, exp, log10
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt


plt.rcParams.update({'font.size': 16})


def frictionCoeff(gravAccel, densityMelt, densityDroplet, diameter, viscosityMelt):
    f = (pi/6) * ((densityDroplet - densityMelt) / densityMelt) * ((densityMelt / viscosityMelt)**2) * gravAccel * (diameter**3)
    return f

def rFromWeberTurbulent(density_melt, density_droplet, gravity, Cd=0.2, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = sqrt(weber_number * (((16 * gravity * (density_droplet - density_melt)) / (3 * Cd * surface_tension)) *
                        ((density_droplet - density_melt) / density_melt))**(-1))
    return r

def rFromWeberLaminar(dynamic_viscosity, density_droplet, density_melt, gravity, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = ((81 * weber_number * surface_tension * (dynamic_viscosity**2)) / (8 * (gravity**2) * ((density_droplet - density_melt)**3)))**(1/5)
    return r

gravity_vesta = 0.25
gravity_earth = 9.8
density_droplet = 7800
density_melt = 3750
diffusion_coeff = 10**-8
vesta_droplet_radius = 0.0185
earth_droplet_radius = 0.002957762

etas = [10 ** (-3.5), 10 ** (-3.0), 10 ** (-2.5), 10 ** (-2.0), 10 ** (-1.5), 10 ** (-1.0)]
radii = [0.001, 0.01, 0.1]

vesta_etas_f = []
vesta_rad_f_1 = []
vesta_rad_f_2 = []
earth_etas_f = []
earth_rad_f_1 = []
earth_rad_f_2 = []


for i in etas:
    vesta_f = frictionCoeff(gravAccel=gravity_vesta, densityMelt=density_melt, densityDroplet=density_droplet,
                            diameter=(2 * vesta_droplet_radius), viscosityMelt=i)
    earth_f = frictionCoeff(gravAccel=gravity_earth, densityMelt=density_melt, densityDroplet=density_droplet,
                            diameter=(2 * earth_droplet_radius), viscosityMelt=i)
    vesta_etas_f.append(vesta_f)
    earth_etas_f.append(earth_f)
for i in radii:
    vesta_f = frictionCoeff(gravAccel=gravity_vesta, densityMelt=density_melt, densityDroplet=density_droplet,
                            diameter=(2 * i), viscosityMelt=(10**-3.5))
    earth_f = frictionCoeff(gravAccel=gravity_earth, densityMelt=density_melt, densityDroplet=density_droplet,
                            diameter=(2 * i), viscosityMelt=(10**-3.5))
    vesta_f_2 = frictionCoeff(gravAccel=gravity_vesta, densityMelt=density_melt, densityDroplet=density_droplet,
                            diameter=(2 * i), viscosityMelt=(10 ** -1.0))
    earth_f_2 = frictionCoeff(gravAccel=gravity_earth, densityMelt=density_melt, densityDroplet=density_droplet,
                            diameter=(2 * i), viscosityMelt=(10 ** -1.0))
    
    vesta_rad_f_1.append(vesta_f)
    vesta_rad_f_2.append(vesta_f_2)
    earth_rad_f_1.append(earth_f)
    earth_rad_f_2.append(earth_f_2)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(etas, [log10(i) for i in vesta_etas_f], linewidth=2.0, color='black', linestyle="-", label="Vesta ($r = 1.85 cm$)")
ax1.plot(etas, [log10(i) for i in earth_etas_f], linewidth=2.0, color='black', linestyle="--", label="Earth ($r = 0.29 cm$)")
ax1.axhline(log10(10), color='black', linestyle="--", label='Turbulent-Laminar Transition')
ax1.set_xscale('log')
ax1.set_xlabel("$\eta (Pa \cdot s)$")
ax1.set_ylabel("log(f)")
ax1.set_title("Friction Coefficient (f) For Stable Radius Metal Droplets on Vesta and Earth")
ax1.grid()
ax1.legend(loc='upper right')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot([i * 100 for i in radii], [log10(i) for i in vesta_rad_f_1], linewidth=2.0, linestyle="-", label="Vesta ($\eta$=10$^{-3.5} Pa \cdot s$)")
ax2.plot([i * 100 for i in radii], [log10(i) for i in vesta_rad_f_2], linewidth=2.0, linestyle="-", label="Vesta ($\eta$=10$^{-1.0} Pa \cdot s$)")
ax2.plot([i * 100 for i in radii], [log10(i) for i in earth_rad_f_1], linewidth=2.0, linestyle="-", label="Earth ($\eta$=10$^{-3.5} Pa \cdot s$)")
ax2.plot([i * 100 for i in radii], [log10(i) for i in earth_rad_f_2], linewidth=2.0, linestyle="-", label="Earth ($\eta$=10$^{-1.0} Pa \cdot s$)")
ax2.axhline(log10(10), color='black', linestyle="--", label='Turbulent-Laminar Transition')
ax2.set_xlabel("Radius (cm)")
ax2.set_ylabel("log(f)")
ax2.set_title("Friction Coefficient (f) For Variable Iron Droplet Radii Over Static $\eta$ on Vesta and Earth")
ax2.grid()
ax2.legend(loc='lower right')






plt.show()