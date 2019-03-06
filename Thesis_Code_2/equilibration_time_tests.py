from math import pi, log, sqrt, log10, exp
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt



def eqTime(A, B, C_m_t, C_m_0):
    t = (1 / A) * log((C_m_t + (B/A)) / (C_m_0 + (B/A)))
    return t

def eqDistance(eqTime, velocity):
    d = velocity * eqTime
    return d

def coeffA(velocity, radius, F_m, F_s, D_ms):
    A = (velocity / (2 * radius)) * ((F_m / (F_m + (F_s / D_ms))) - 1)
    return A

def coeffB(velocity, radius, F_m, F_s, D_ms, C_0_s):
    B = (velocity / (2 * radius)) * ((F_s * C_0_s) / (F_m + (F_s / D_ms)))
    return B

def calcFm(radius, density_droplet, density_melt, D_s, velocity):
    vol_metal = (radius**3) * density_droplet / 3
    vol_silicate = (radius**2) * density_melt * sqrt((2 * D_s * radius) / velocity)
    Fm = vol_metal / (vol_metal + vol_silicate)
    return Fm

def calcFs(Fm):
    Fs = 1 - Fm
    return Fs

def surfaceGravity(radius, melt_density):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius**3 * melt_density) / radius**2
    return gSurface

def gravity(radius, melt_density, depth):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius ** 3 * melt_density) / radius**2
    gAccel = gSurface * (1 - (depth/radius))
    return gAccel

def adiabat(current_depth, current_temperature, depth_increment, depth, gravity, thermal_expansion, heat_capacity,
            depth_gradient, adiabat_gradient):
    current_depth += depth_increment
    # dT/dh = (alpha * g * T) / c_p
    current_temperature += ((thermal_expansion * gravity * current_temperature) / heat_capacity) * depth_increment
    adiabat_gradient.append(current_temperature)
    depth_gradient.append(current_depth / 1000)
    if current_depth < depth:  # recursion!
        return adiabat(current_depth, current_temperature, depth_increment, depth, gravity, thermal_expansion,
                       heat_capacity, depth_gradient, adiabat_gradient)
    else:
        return depth_gradient, adiabat_gradient

def hydrostatic(current_depth, depth_increment, depth, gravity, density_melt, current_pressure, depth_gradient,
                pressure_gradient):
    current_depth += depth_increment
    # dP/dh = rho * g
    current_pressure += (density_melt * gravity) * depth_increment
    pressure_gradient.append(current_pressure * (10**(-9)))  # gPa
    depth_gradient.append(current_depth / 1000)
    if current_depth < depth:
        return hydrostatic(current_depth, depth_increment, depth, gravity, density_melt, current_pressure, depth_gradient,
                pressure_gradient)
    else:
        return depth_gradient, pressure_gradient

def turbulentVelocity(gravAccel, radius_droplet, density_droplet, density_melt):
    Cd = 0.2
    v = sqrt(((4 * gravAccel * (2 * radius_droplet)) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v

def laminarVelocity(density_droplet, density_melt, gravity, droplet_radius, dynamic_viscosity):
    v = ((density_droplet - density_melt) * gravity * ((2 * droplet_radius)**2)) / (18 * dynamic_viscosity)
    return v


def rubieCmt(A, B, time, C_m_0):
    if C_m_0 > 0:
        c_m_t = ((C_m_0 + (B / A)) * exp(A * time)) - (B / A)
        return c_m_t
    else:
        return 0



if __name__ == "__main__":

    radius_vesta = 265 * 1000  # 265 km
    current_vesta_core_radius_min = 107 * 1000
    current_vesta_core_radius_max = 113 * 1000
    depth_increment = 5 * 1000  # 5 km
    droplet_radius = 0.01  # m
    densityDroplet = 7800  # kg m3
    densityMelt = 3750  # kg m3
    c_0_m = 32700  # ppm
    c_0_s = 42  # ppm
    chemical_diffusivity = 10 ** (-8)
    thermal_expansion = 6 * 10 ** (-5)
    heat_capacity = 10 ** 3
    dynamic_viscosity = 10**(-2)
    thermal_diffusivity = 10 ** (-6)
    distrib_coeff = 28
    surface_gravity = surfaceGravity(radius=radius_vesta, melt_density=densityMelt)

    adiabatic_depths, adiabatic = adiabat(
        current_depth=0 / 1000,  # begin at surface
        current_temperature=2000,  # degrees K
        depth_increment=depth_increment,  # 10 km interval
        depth=radius_vesta,  # to depth of 250 km
        gravity=surface_gravity,  # m/s^2, Vesta
        thermal_expansion=6 * 10 ** (-5),
        heat_capacity=10 ** 3,
        depth_gradient=[0],  # returns this list at index=0
        adiabat_gradient=[2000],  # returns this list at index=1
    )

    hydrostatic_depths, hydrostat = hydrostatic(
        current_depth=0 / 1000,
        depth_increment=depth_increment,
        depth=radius_vesta,
        gravity=surface_gravity,
        current_pressure=0,
        density_melt=densityMelt,
        depth_gradient=[0],
        pressure_gradient=[0],
    )

    F_m_list = []
    F_s_list = []
    A_list = []
    B_list = []
    eq_time_list = []
    c_m_t_dict = {}

    rads = np.arange(0.001, 0.015, 0.002)
    times = np.arange(0, 1000, 1)

    for r in rads:
        velocity = turbulentVelocity(gravAccel=surface_gravity, density_droplet=densityDroplet, density_melt=densityMelt, radius_droplet=r)
        F_m = calcFm(radius=r, D_s=chemical_diffusivity, density_droplet=densityDroplet, density_melt=densityMelt, velocity=velocity)
        F_s = calcFs(Fm=F_m)
        A = coeffA(velocity=velocity, D_ms=distrib_coeff, F_m=F_m, F_s=F_s, radius=r)
        B = coeffB(velocity=velocity, D_ms=distrib_coeff, F_m=F_m, F_s=F_s, radius=r, C_0_s=c_0_s)
        F_m_list.append(F_m)
        F_s_list.append(F_s)


        A_list.append(A)
        B_list.append(B)

        c_m_t_list = []
        C_0_s_t = c_0_s
        for t in times:
            B_adj = coeffB(velocity=velocity, D_ms=distrib_coeff, F_m=F_m, F_s=F_s, radius=r, C_0_s=C_0_s_t)
            c_m_t = rubieCmt(A=A, B=B, C_m_0=C_0_s_t, time=t)
            c_m_t_list.append(c_m_t)
        c_m_t_dict.update({str(r): c_m_t_list})

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot([i * 100 for i in rads], F_m_list, linewidth=2, label="Mass Fraction Metal Equilibrated")
    ax1_2 = ax1.twinx()
    ax1_2.plot(rads, F_s_list, linewidth=2, linestyle="--", label="Mass Fraction Silicate Equilibrated")
    ax1.grid()
    ax1.set_xlabel("Droplet Radius (cm)")
    ax1.set_ylabel("Mass Fraction Equilibrated Metal (Solid)")
    ax1_2.set_ylabel("Mass Fraction Equilibrated Silicate Melt (Dashed)")
    ax1.set_title("Mass Fraction Equilibrated Metal vs Mass Fraction Equilibrated Silicate Melt (Fm + Fs = 1)")


    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(A_list, [i * 100 for i in rads], linewidth=2, label="Coefficient A vs Droplet Radius (solid)")
    ax2_2 = ax2.twinx()
    ax2_2.plot(A_list, F_m_list, linewidth=2, linestyle="--", label="Coefficient A vs Fm (dashed)")
    ax2.grid()
    ax2.set_xlabel("Variable A")
    ax2.set_ylabel("Droplet Radius (cm) (solid)")
    ax2_2.set_ylabel("Equilibrated Metal Mass Fraction (Fm) (dashed)")
    ax2.set_title("Variable A vs Fm and Droplet Radius (Rubie et al 2003)")

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(B_list, [i * 100 for i in rads], linewidth=2, label="Coefficient B vs Droplet Radius (solid)")
    ax3_2 = ax2.twinx()
    ax3_2.plot(B_list, F_m_list, linewidth=2, linestyle="--", label="Coefficient B vs Fm (dashed)")
    ax3.grid()
    ax3.set_xlabel("Variable B")
    ax3.set_ylabel("Droplet Radius (cm) (solid)")
    ax3_2.set_ylabel("Equilibrated Metal Mass Fraction (Fm) (dashed)")
    ax3.set_title("Variable B vs Fm and Droplet Radius (Rubie et al 2003)")

    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111)
    ax4.plot(A_list, B_list, linewidth=2)
    ax4.grid()
    ax4.set_xlabel("Coefficient A")
    ax4.set_ylabel("Coefficient B")
    ax4.set_title("Coefficient A vs Coefficient B")

    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111)
    ax5.plot([i * 100 for i in rads], A_list, linewidth=2, label="Variable A")
    ax5_2 = ax5.twinx()
    ax5_2.plot([i * 100 for i in rads], B_list, linewidth=2, linestyle="--", label="Variable B")
    ax5.grid()
    ax5.set_xlabel("Droplet Radius (cm)")
    ax5.set_ylabel("Variable A (solid)")
    ax5_2.set_ylabel("Variable B (dashed)")
    ax5.set_title("Droplet Radius Variables A and B")


    fig6 = plt.figure()
    ax6 = fig6.add_subplot(111)
    ax6.plot([i * 100 for i in rads], [x/y for x, y in zip(B_list, A_list)], linewidth=2)
    # ax6.plot([i * 100 for i in rads], [round(x / y, 6) for x, y in zip(B_list, A_list)], linewidth=2)
    ax6.set_xlabel("Droplet Radius (cm)")
    ax6.set_ylabel("B/A")
    ax6.set_title("Droplet Radius vs. B/A (B/A = minimum droplet concentration) (Rubie et al 2003)")
    ax6.grid()

    fig7 = plt.figure()
    ax7 = fig7.add_subplot(111)
    ax7.plot(F_m_list, [x/y for x, y in zip(B_list, A_list)], linewidth=2)
    # ax7.plot(F_m_list, [round(x / y, 6) for x, y in zip(B_list, A_list)], linewidth=2)
    ax7.set_xlabel("Mass Fraction Metal Equilibrated (Fm)")
    ax7.set_ylabel("B/A")
    ax7.set_title("Mass Fraction Metal Equilibrated vs. B/A (B/A = minimum droplet concentration) (Rubie et al 2003)")
    ax7.grid()

    fig8 = plt.figure()
    ax8 = fig8.add_subplot(111)
    for key in c_m_t_dict.keys():
        ax8.plot(times, c_m_t_dict[key], linewidth=2, label="Droplet Radius: {} cm".format(round(float(key) * 100, 4)))
    ax8.set_xlabel("Time (s)")
    ax8.set_ylabel("Concentration of Droplet (C_m_t)")
    ax8.set_title("Droplet Concentration Over Time (Rubie et al 2003)")
    ax8.grid()
    ax8.legend(loc='lower right')

    print([x/y for x, y in zip(B_list, A_list)])

    plt.show()
