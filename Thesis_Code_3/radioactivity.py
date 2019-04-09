import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rcParams.update({'font.size': 16})


def decay(half_life, curr_nuclii, max_time, timestep, current_time, original_nuclii, rad_list=[]):
    decay_const = log(0.5) / half_life
    if current_time <= max_time:
        remaining_nuclii = curr_nuclii * exp(decay_const * timestep)
        rad_list.append((remaining_nuclii / original_nuclii))
        return decay(half_life=half_life, curr_nuclii=remaining_nuclii, max_time=max_time, timestep=timestep,
                     current_time=current_time + timestep, original_nuclii=original_nuclii, rad_list=rad_list)
    else:
        return rad_list

def avg_vals_time(list_of_lists):
    avgs = []
    z_list = list(zip(*list_of_lists))
    for i in z_list:
        avg = sum(i) / len(i)
        avgs.append(avg)
    return avgs


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
My_5_index = time_list.index(5)

hf_decay = decay(half_life=hf_half_life, curr_nuclii=original_hf, max_time=max_time, timestep=timestep,
                 current_time=timestep, original_nuclii=original_hf, rad_list=[original_hf / original_hf])
al_decay = decay(half_life=al_half_life, curr_nuclii=original_al, max_time=max_time, timestep=timestep,
                 current_time=timestep, original_nuclii=original_al,  rad_list=[original_al / original_al])
fe_decay = decay(half_life=fe_half_life, curr_nuclii=original_fe, max_time=max_time, timestep=timestep,
                 current_time=timestep, original_nuclii=original_fe,  rad_list=[original_fe / original_fe])

w_abundance = [1 - i for i in hf_decay]

hf_rel_at_5My = hf_decay[My_5_index]
al_rel_at_5My = al_decay[My_5_index]
fe_rel_at_5My = fe_decay[My_5_index]


print(1 - hf_rel_at_5My, 1 - al_rel_at_5My, 1 - fe_rel_at_5My)
print(hf_rel_at_5My, al_rel_at_5My / original_al, fe_rel_at_5My / original_fe)


hf_at_5 = hf_decay[My_5_index]
al_at_5 = al_decay[My_5_index]
fe_at_5 = fe_decay[My_5_index]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(time_list, hf_decay, linewidth=2.0, label='182Hf (hl = 8.9 * 10^6 y)')
ax.plot(time_list, al_decay, linewidth=2.0, label='26Al (hl = 7.17 * 10^5 y)')
ax.plot(time_list, fe_decay, linewidth=2.0, label='60Fe (hl = 2.6 * 10^6 y)')
ax.plot(time_list, w_abundance, linewidth=2.0, linestyle="--", label="182W (stable)")
ax.axvspan(0, 5, alpha=0.2, color='red')
ax.set_xlabel("Time (My)")
ax.set_ylabel("Relative Isotope Abundance (remaining/original)")
ax.set_title("Isotope Decay")
ax.legend(loc='center right')
minorLocator = MultipleLocator((timestep * 5) / 10**6)
ax.xaxis.set_minor_locator(minorLocator)
ax.grid()

eucrite_df = pd.read_excel("eucrites_kleine_2002.xlsx")

sample_name_list = []
w_182_w_184_list = []  # absolute ratio over time
w_182_w_184_eucrite_relative_list = []  # relative ratio to the final measured eucrite abundance
hf_180_w_184_list = []
epsilon_w_list = []

My_5_w_182_w_184 = w_abundance[My_5_index]

for row in eucrite_df.index:
   sample_name = eucrite_df['Sample'][row]
   w_182_w_184 = eucrite_df['182W/184W'][row]
   hf_180_w_184 = eucrite_df['180Hf/184W'][row]
   epsilon_w = eucrite_df['epsilon_W'][row]

   w_182_w_184_time = [i * float(w_182_w_184) for i in w_abundance]
   w_182_w_184_time_rel = [i / float(w_182_w_184) for i in w_182_w_184_time]
   epsilon_w_time = [((i / w_182_w_184_terrestrial) - 1) * (10**4) for i in w_182_w_184_time]

   sample_name_list.append(sample_name)
   w_182_w_184_list.append(w_182_w_184_time)
   w_182_w_184_eucrite_relative_list.append(My_5_w_182_w_184)
   epsilon_w_list.append(epsilon_w_time)

w_182_w_184_list_avgs = avg_vals_time(list_of_lists=w_182_w_184_list)


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.axvspan(0, 5, alpha=0.2, color='red')
ax2.set_xlabel("Time (My)")
ax2.set_ylabel("Relative 182W/184W")
ax2.set_title("182W/184W In Eucrites Over Time")
minorLocator = MultipleLocator((timestep * 5) / 10**6)
ax2.xaxis.set_minor_locator(minorLocator)
ax2.grid()
for index, i in enumerate(sample_name_list):
    ax2.plot(time_list, w_182_w_184_list[index], label="182W/184W ({})".format(i))
    ax2.axhline(w_182_w_184_list_avgs[My_5_index], linestyle="--", color='black')
    ax2.annotate("Avg. 182W/184W (5 My) = {}".format(round(float(w_182_w_184_list_avgs[My_5_index]), 6)),
                 (time_list[My_5_index], w_182_w_184_list_avgs[My_5_index]), xytext=(20.2, 0.42),
                 arrowprops=dict(facecolor='black', shrink=0.05))
    ax2.annotate("Avg. 182W/184W (100 My) = {}".format(round(float(w_182_w_184_list_avgs[-1]), 6)),
                 (time_list[-1], w_182_w_184_list_avgs[-1]), xytext=(70.2, 0.62),
                 arrowprops=dict(facecolor='black', shrink=0.05))

ax2.plot(time_list, w_182_w_184_list_avgs, label="Average 182W/184W", linestyle="--")

ax2.legend(loc='lower right')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.axvspan(0, 5, alpha=0.2, color='red')
ax3.set_xlabel("Time (My)")
ax3.set_ylabel("Epsilon 182W")
ax3.set_title("Epsilon 182W In Eucrites Over Time")
minorLocator = MultipleLocator((timestep * 5) / 10**6)
ax3.xaxis.set_minor_locator(minorLocator)
ax3.grid()

for index, i in enumerate(sample_name_list):
    ax3.plot(time_list, epsilon_w_list[index], label="Epsilon 182W ({})".format(i))

ax3.legend(loc='center right')

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.axvspan(0, 5, alpha=0.2, color='red')
ax4.set_xlabel("Time (My)")
ax4.set_ylabel("Relative $^{182}$W/$^{184}$W")
ax4.set_title("Relative $^{182}$W/$^{184}$W In Eucrites Over Time")
minorLocator = MultipleLocator((timestep * 5) / 10**6)
ax4.xaxis.set_minor_locator(minorLocator)
ax4.grid()
for index, i in enumerate(sample_name_list):
    ax4.plot(time_list, [j / w_182_w_184_list[index][-1] for j in w_182_w_184_list[index]], label="182W/184W ({})".format(i))
    ax4.axhline(w_182_w_184_list_avgs[My_5_index] / w_182_w_184_list[index][-1], linestyle="--", color='black')
    ax4.annotate("Avg. 182W/184W (5 My) = {}".format(round(float(w_182_w_184_list_avgs[My_5_index] / w_182_w_184_list_avgs[-1]), 6)),
                 (time_list[My_5_index], w_182_w_184_list_avgs[My_5_index] / w_182_w_184_list_avgs[-1]), xytext=(20.2, 0.42),
                 arrowprops=dict(facecolor='black', shrink=0.05))
    ax4.annotate("Avg. 182W/184W (100 My) = {}".format(round(float(w_182_w_184_list_avgs[-1] / w_182_w_184_list_avgs[-1]), 6)),
                 (time_list[-1], w_182_w_184_list_avgs[-1] / w_182_w_184_list_avgs[-1]), xytext=(80.2, 0.82),
                 arrowprops=dict(facecolor='black', shrink=0.05))

ax4.plot(time_list, [j / w_182_w_184_list_avgs[-1] for j in w_182_w_184_list_avgs], label="Average 182W/184W", linestyle="--")

ax4.legend(loc='center right')

print("***{}".format(w_182_w_184_list_avgs[My_5_index]))

# fig4 = plt.figure()
# ax4 = fig4.add_subplot(111)
# ax4.axvspan(0, 5, alpha=0.2, color='red')
# ax4.set_xlabel("Time (My)")
# ax4.set_ylabel("Relative 182W/184W")
# ax4.set_title("182W/184W In Eucrites Over Time")
# minorLocator = MultipleLocator((timestep * 5) / 10**6)
# ax4.xaxis.set_minor_locator(minorLocator)
# ax4.grid()
#
# for index, i in enumerate(sample_name_list):
#     ax4.plot(time_list, w_182_w_184_eucrite_relative_list[index], label="182W/184W ({})".format(i))
#
# ax4.legend(loc='center right')

plt.show()

