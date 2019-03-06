# import matplotlib as mpl
# mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def decay(half_life, curr_nuclii, max_time, timestep, current_time, original_nuclii, rad_list=[]):
    decay_const = log(0.5) / half_life
    if current_time <= max_time:
        remaining_nuclii = curr_nuclii * exp(decay_const * timestep)
        rad_list.append((remaining_nuclii / original_nuclii))
        return decay(half_life=half_life, curr_nuclii=remaining_nuclii, max_time=max_time, timestep=timestep,
                     current_time=current_time + timestep, original_nuclii=original_nuclii, rad_list=rad_list)
    else:
        return rad_list


hf_half_life = 8.9 * 10**6
al_half_life = 7.17 * 10**5
fe_half_life = 3.0 * 10**5
max_time = 100 * 10**6
original_hf = 100
original_al = 100
original_fe = 100
timestep = 1 * 10**6
time_list = [i / (1 * 10**6) for i in np.arange(0, max_time + timestep, timestep)]

hf_decay = decay(half_life=hf_half_life, curr_nuclii=original_hf, max_time=max_time, timestep=timestep,
                 current_time=timestep, original_nuclii=original_hf, rad_list=[original_hf / original_hf])
al_decay = decay(half_life=al_half_life, curr_nuclii=original_al, max_time=max_time, timestep=timestep,
                 current_time=timestep, original_nuclii=original_al,  rad_list=[original_al / original_al])
fe_decay = decay(half_life=fe_half_life, curr_nuclii=original_fe, max_time=max_time, timestep=timestep,
                 current_time=timestep, original_nuclii=original_fe,  rad_list=[original_fe / original_fe])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(time_list, hf_decay, linewidth=2.0, label='182Hf (hl = 8.9 * 10^6 y)')
ax.plot(time_list, al_decay, linewidth=2.0, label='26Al (hl = 7.17 * 10^5 y)')
ax.plot(time_list, fe_decay, linewidth=2.0, label='60Fe (hl = 3.0 * 10^5 y)')
ax.set_xlabel("Time (my)")
ax.set_ylabel("Relative Isotope Abundance (remaining/original)")
ax.set_title("Isotope Decay")
ax.legend(loc='upper right')
minorLocator = MultipleLocator((timestep * 5) / 10**6)
ax.xaxis.set_minor_locator(minorLocator)
ax.grid()

plt.show()
