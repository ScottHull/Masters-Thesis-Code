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
w_182_w_184_terrestrial = 0.864900  # Kleine & Walker 2017 Tungsten Isotopes in Planets
w_182_w_184_terrestrial_old = 0.864680  # Kleine et al. 2002 Eucrites
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

w_abundance = [1 - i for i in hf_decay]


my_5_index = time_list.index(5)
hf_at_5 = hf_decay[my_5_index]
al_at_5 = al_decay[my_5_index]
fe_at_5 = fe_decay[my_5_index]

eucrite_df = pd.read_excel("eucrites_kleine_2002.xlsx")

sample_name_list = []
w_182_w_184_list = []
hf_180_w_184_list = []
epsilon_w_list = []

for row in eucrite_df.index:
   sample_name = eucrite_df['Sample'][row]
   w_182_w_184 = eucrite_df['182W/184W'][row]
   hf_180_w_184 = eucrite_df['180Hf/184W'][row]
   epsilon_w = eucrite_df['epsilon_W'][row]

   w_182_w_184_time = [i * float(w_182_w_184) for i in w_abundance]
   epsilon_w_time = [((i / w_182_w_184_terrestrial_old) - 1) * (10**4) for i in w_182_w_184_time]

   sample_name_list.append(sample_name)
   w_182_w_184_list.append(w_182_w_184_time)
   epsilon_w_list.append(epsilon_w_time)
