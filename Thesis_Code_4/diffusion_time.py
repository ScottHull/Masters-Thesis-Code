import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, sqrt, pi, log10
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rcParams.update({'font.size': 16})


def diffusionTime(diffusivity, length):
    t = (length**2) / diffusivity

    return t

vesta_radius = 263 * 1000  # m
vesta_core_radius = 113 * 1000  # m
vesta_mantle_depth = vesta_radius - vesta_core_radius
w_diffusivity = 10**(-8)
lengths = list(np.arange(0, 5000, 10))
times_1 = [diffusionTime(diffusivity=w_diffusivity, length=int(i)) for i in lengths]
# times_2 = [diffusionTime(diffusivity=w_diffusivity * (10**-8), length=int(i)) for i in lengths]
times_1_millions_of_years = [(i * (3.17098 * (10**-8))) / (10**6) for i in times_1]
# times_2_millions_of_years = [(i * (3.17098 * (10**-8))) / (10**6) for i in times_2]


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot([i / 1000 for i in lengths], times_1_millions_of_years, linewidth=2.0, color='black', label="$D=10^{-8}$")
# ax.plot([i / 1000 for i in lengths], times_2_millions_of_years, linewidth=2.0, color='green', label="$D=10^{-16}$")
ax.grid()
ax.set_title("W Diffusion Time as a Function of Length")
ax.set_xlabel("Length (km)")
ax.set_ylabel("Time (Ma)")
# ax.axhline(5, linewidth=2.0, linestyle="--", color="red", label="End of Vesta Metal-Silicate Separation")
# ax.axhline(100, linewidth=2.0, linestyle="-.", color="red", label="Complete Decay of $^{182}$Hf")
# ax.axhline(4500, linewidth=2.0, linestyle="--", color='green', label="Present Day")
# ax.axhline(vesta_mantle_depth / 1000, linewidth=2.0, linestyle="--", color='blue', label="Vesta Mantle Depth")
# ax.axvline(vesta_radius / 1000, linewidth=2.0, linestyle="-.", color='blue', label="Vesta Radius")
ax.legend(loc='upper left')


plt.show()
