from math import pi, log, sqrt, exp
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

def frac_eq(radius, diffusivity, settling_velocity, dynamic_viscosity, z, density_melt):
    kinematic_viscosity = dynamic_viscosity / density_melt
    num = -0.996 * (diffusivity**(2/3)) * z
    den = ((2 * radius)**(3/2)) * (settling_velocity**(1/2)) * (kinematic_viscosity**(1/6))
    f = 1 - exp(num / den)
    return f


etas = [10**-3.5, 10**-1.0]
# radii = list(np.arange(0.001, 0.02 + 0.001, 0.001))
radii = [0.001, 0.01, 0.1, 1, 10, 100, 1000]
gravity_vesta = 0.25
gravity_earth = 9.8
density_droplet = 7800
density_melt = 3750
diffusion_coeff = 10**-8
earth_magma_ocean_depth = 500 * 1000
radius_vesta = 263 * 1000
radius_vesta_core = 113 * 1000
radius_vesta_mantle = radius_vesta - radius_vesta_core

velocities_vesta = []
z_eqs_vesta = []
velocities_earth = []
z_eqs_earth = []
frac_eq_vesta_list = []
frac_eq_earth_list = []

for eta in etas:
    velocities_vesta_temp = []
    z_eqs_vesta_temp = []
    velocities_earth_temp = []
    z_eqs_earth_temp = []
    frac_eq_vesta_list_temp = []
    frac_eq_earth_list_temp = []
    for i in radii:
        v_earth = turbulentVelocity(gravity=gravity_earth, density_melt=density_melt, density_droplet=density_droplet,
                                    droplet_radius=i)
        z_earth = z_eq(radius_droplet=i, dynamic_viscosity=eta, settling_velocity=v_earth, diffusion_coeff=diffusion_coeff,
                       density_melt=density_melt)
        v_vesta = turbulentVelocity(gravity=gravity_vesta, density_melt=density_melt, density_droplet=density_droplet,
                                    droplet_radius=i)
        z_vesta = z_eq(radius_droplet=i, dynamic_viscosity=eta, settling_velocity=v_vesta, diffusion_coeff=diffusion_coeff,
                       density_melt=density_melt)
        frac_eq_vesta = frac_eq(radius=i, dynamic_viscosity=eta, settling_velocity=v_vesta, diffusivity=diffusion_coeff,
                                z=radius_vesta_mantle, density_melt=density_melt)
        frac_eq_earth = frac_eq(radius=i, dynamic_viscosity=eta, settling_velocity=v_earth, diffusivity=diffusion_coeff,
                                z=earth_magma_ocean_depth, density_melt=density_melt)
        velocities_vesta_temp.append(v_vesta)
        z_eqs_vesta_temp.append(z_vesta)
        velocities_earth_temp.append(v_earth)
        z_eqs_earth_temp.append(z_earth)
        frac_eq_vesta_list_temp.append(frac_eq_vesta)
        frac_eq_earth_list_temp.append(frac_eq_earth)
    velocities_vesta.append(velocities_vesta_temp)
    z_eqs_vesta.append(z_eqs_vesta_temp)
    velocities_earth.append(velocities_earth_temp)
    z_eqs_earth.append(z_eqs_earth_temp)
    frac_eq_vesta_list.append(frac_eq_vesta_list_temp)
    frac_eq_earth_list.append(frac_eq_earth_list_temp)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot([i for i in radii], z_eqs_vesta[0], linewidth=2.0, color='blue', label="Vesta (g=0.25 m/s$^2$, $\eta$=10$^{-3.5}$)")
ax1.plot([i for i in radii], z_eqs_earth[0], linewidth=2.0, color='red', label="Earth (g=9.80 m/s$^2$, $\eta$=10$^{-3.5}$)")
ax1.plot([i for i in radii], z_eqs_vesta[1], linewidth=2.0, color='blue', linestyle="--", label="Vesta (g=0.25 m/s$^2$, $\eta$=10$^{-1.0}$)")
ax1.plot([i for i in radii], z_eqs_earth[1], linewidth=2.0, color='red', linestyle="--", label="Earth (g=9.80 m/s$^2$, $\eta$=10$^{-1.0}$)")
ax1.axhline(2891 * 1000, linestyle="--", color='black', label="Earth Mantle Depth")
ax1.axhline(149700, linestyle="-.", color='black', label="Vesta Mantle Depth")
ax1.set_xlabel("Droplet Radius (m)")
ax1.set_ylabel("z$_{eq}$ (m)")
ax1.set_title("Droplet Radius vs. 99% Chemical Equilibration Distance")
ax1.grid()
ax1.legend(loc='upper left')
ax1.set_yscale('log')
ax1.set_xscale('log')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot([i for i in radii], frac_eq_vesta_list[0], linewidth=2.0, color='blue', label="Vesta (g=0.25 m/s$^2$, $\eta$=10$^{-3.5}$)")
ax2.plot([i for i in radii], frac_eq_earth_list[0], linewidth=2.0, color='red', label="Earth (g=9.80 m/s$^2$, $\eta$=10$^{-3.5}$)")
ax2.plot([i for i in radii], frac_eq_vesta_list[1], linewidth=2.0, color='blue', linestyle="--", label="Vesta (g=0.25 m/s$^2$, $\eta$=10$^{-1.0}$)")
ax2.plot([i for i in radii], frac_eq_earth_list[1], linewidth=2.0, color='red', linestyle="--", label="Earth (g=9.80 m/s$^2$, $\eta$=10$^{-1.0}$)")
ax2.set_xlabel("Droplet Radius (m)")
ax2.set_ylabel("Fraction Equilibrated (%)")
ax2.set_title("Fraction of Metal-Silicate Chemical Equilibrium as a Function of Droplet Radius")
ax2.grid()
ax2.legend(loc='lower left')
ax2.set_xscale('log')
ax2.set_yscale('log')

plt.show()

