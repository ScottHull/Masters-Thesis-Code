import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp
import numpy as np


def decayModel(half_life, curr_nuclii, max_time, timestep, current_time=0, N_t_list=[]):
    decay_const = log(0.5) / half_life
    if current_time <= max_time:
        remaining_nuclii = curr_nuclii * exp(-decay_const * timestep)
        N_t_list.append(remaining_nuclii)
        return decayModel(half_life=half_life, curr_nuclii=remaining_nuclii, max_time=max_time, timestep=timestep,
                          current_time=current_time+timestep, N_t_list=N_t_list)
    else:
        return N_t_list


hf_182_half_life = 8.9 * 10**6
original_num_nucleii = 100

max_time = 100 * 10**6
timestep = 1 * 10**6
timesteps = np.arange(0, max_time + timestep, timestep)

hf_decay = decayModel(half_life=hf_182_half_life, curr_nuclii=original_num_nucleii, max_time=max_time, timestep=timestep)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(timesteps, hf_decay)
ax.grid()



plt.show()
