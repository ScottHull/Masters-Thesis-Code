import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes
import pandas as pd
from math import log, exp, sqrt, pi, log10
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rcParams.update({'font.size': 16})




model_epsilon_w182_5ma_df = pd.read_csv("epsilon_182w_5ma.csv")
model_epsilon_w182_100ma_df = pd.read_csv("epsilon_182w_100ma.csv")
data_epsilon_w182_df = pd.read_excel("epsilon_182w_values.xlsx").set_index("Sample")
# fO2 = [3.5, 2.0, 1.0, 0.5, 0, -0.8, -1.10, -2.25, -2.45, -3.5]
fO2 = [-0.8, -1.10, -2.25, -2.45]
components = ['core', 'mantle', 'bulk']

iiab_iron_meteorite = data_epsilon_w182_df['new epsilon 182']['IIAB Iron Meteorite']
ivb_iron_meteorite = data_epsilon_w182_df['new epsilon 182']['IVB Iron Meteorite']
chondrite = data_epsilon_w182_df['new epsilon 182']['Chondrite']
solar_system_initial = data_epsilon_w182_df['new epsilon 182']['Solar System Initial']
eucrite = data_epsilon_w182_df['new epsilon 182']['Eucrite']

epsilon_w182_5ma_dict = {}
epsilon_w182_100ma_dict = {}

for i in components:
    epsilon_w182_5ma_dict.update({i: {'fO2': [], 'epsilon_182w': [], 'name': []}})
    epsilon_w182_100ma_dict.update({i: {'fO2': [], 'epsilon_182w': [], 'name': []}})

for i in components:
    for j in fO2:
        header = "{}-fO2_{}".format(i, j)
        epsilon_w182_5ma = model_epsilon_w182_5ma_df[header][0]
        epsilon_w182_5ma_dict[i]['fO2'].append(j)
        epsilon_w182_5ma_dict[i]['epsilon_182w'].append(epsilon_w182_5ma)
        epsilon_w182_5ma_dict[i]['name'].append('$\Delta$IW' + str(j))
        epsilon_w182_100ma = model_epsilon_w182_100ma_df[header][0]
        epsilon_w182_100ma_dict[i]['fO2'].append(j)
        epsilon_w182_100ma_dict[i]['epsilon_182w'].append(epsilon_w182_100ma)
        epsilon_w182_100ma_dict[i]['name'].append('$\Delta$IW' + str(j))

epsilon_w182_5ma_dict['core']['epsilon_182w'].append(iiab_iron_meteorite)
epsilon_w182_5ma_dict['core']['name'].append("IIAB")
epsilon_w182_5ma_dict['core']['epsilon_182w'].append(ivb_iron_meteorite)
epsilon_w182_5ma_dict['core']['name'].append("IVB")
epsilon_w182_5ma_dict['mantle']['epsilon_182w'].append(eucrite)
epsilon_w182_5ma_dict['mantle']['name'].append("Eucrite")
epsilon_w182_5ma_dict['bulk']['epsilon_182w'].append(solar_system_initial)
epsilon_w182_5ma_dict['bulk']['name'].append("SSI")
epsilon_w182_5ma_dict['bulk']['epsilon_182w'].append(chondrite)
epsilon_w182_5ma_dict['bulk']['name'].append("Chondrite")

epsilon_w182_100ma_dict['core']['epsilon_182w'].append(iiab_iron_meteorite)
epsilon_w182_100ma_dict['core']['name'].append("IIAB")
epsilon_w182_100ma_dict['core']['epsilon_182w'].append(ivb_iron_meteorite)
epsilon_w182_100ma_dict['core']['name'].append("IVB")
epsilon_w182_100ma_dict['mantle']['epsilon_182w'].append(eucrite)
epsilon_w182_100ma_dict['mantle']['name'].append("Eucrite")
epsilon_w182_100ma_dict['bulk']['epsilon_182w'].append(solar_system_initial)
epsilon_w182_100ma_dict['bulk']['name'].append("SSI")
epsilon_w182_100ma_dict['bulk']['epsilon_182w'].append(chondrite)
epsilon_w182_100ma_dict['bulk']['name'].append("Chondrite")

# y_val = 2
#
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# ax1.yaxis.set_major_formatter(plt.NullFormatter())
# ax1.yaxis.set_ticks_position('none')
# ax1.xaxis.grid(True)
# ax1.set_title("$\epsilon$$^{182}$W at Time $t=5 Ma$")
# ax1.set_xlabel("$\epsilon$$^{182}$W")
# for index, i in enumerate(epsilon_w182_5ma_dict['core']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['core']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['core']['name'][index]
#     ax1.scatter(epsilon, y_val, marker='^', s=80, color='blue')
#     ax1.annotate(name, (epsilon + 120, y_val - 0.4))
#     y_val += 2
# for index, i in enumerate(epsilon_w182_5ma_dict['mantle']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['mantle']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['mantle']['name'][index]
#     ax1.scatter(epsilon, y_val, marker='s', s=80, color='green')
#     ax1.annotate(name, (epsilon + 120, y_val - 0.4))
#     y_val += 2
# for index, i in enumerate(epsilon_w182_5ma_dict['bulk']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['bulk']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['bulk']['name'][index]
#     ax1.scatter(epsilon, y_val, marker='o', s=80, color='red')
#     ax1.annotate(name, (epsilon + 120, y_val - 0.4))
#     y_val += 2

# y_val = 2

# f,(ax1,ax2) = plt.subplots(1, 2, sharey=True, facecolor='w')
# # ax1.yaxis.set_major_formatter(plt.NullFormatter())
# # ax1.yaxis.set_ticks_position('none')
# # ax1.xaxis.grid(True)
# ax1.set_title("$\epsilon$$^{182}$W at Time $t=5 Ma$")
# ax1.set_xlabel("$\epsilon$$^{182}$W")
# for index, i in enumerate(epsilon_w182_5ma_dict['core']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['core']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['core']['name'][index]
#     ax1.scatter(epsilon, y_val, marker='^', s=80, color='blue')
#     ax1.annotate(name, (epsilon + 120, y_val - 0.4))
#     ax2.scatter(epsilon, y_val, marker='^', s=80, color='blue')
#     ax2.annotate(name, (epsilon + 120, y_val - 0.4))
#     y_val += 2
# for index, i in enumerate(epsilon_w182_5ma_dict['mantle']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['mantle']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['mantle']['name'][index]
#     ax1.scatter(epsilon, y_val, marker='s', s=80, color='green')
#     ax1.annotate(name, (epsilon + 120, y_val - 0.4))
#     ax2.scatter(epsilon, y_val, marker='s', s=80, color='green')
#     ax2.annotate(name, (epsilon + 120, y_val - 0.4))
#     y_val += 2
# for index, i in enumerate(epsilon_w182_5ma_dict['bulk']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['bulk']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['bulk']['name'][index]
#     ax1.scatter(epsilon, y_val, marker='o', s=80, color='red')
#     ax1.annotate(name, (epsilon + 120, y_val - 0.4))
#     ax2.scatter(epsilon, y_val, marker='o', s=80, color='red')
#     ax2.annotate(name, (epsilon + 120, y_val - 0.4))
#     y_val += 2
#
# ax1.spines['right'].set_visible(False)
# ax2.spines['left'].set_visible(False)
# ax1.yaxis.tick_left()
# ax1.tick_params(labelright='off')
# ax2.yaxis.tick_right()
# d = .015 # how big to make the diagonal lines in axes coordinates
# # arguments to pass plot, just so we don't keep repeating them
# kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
# ax1.plot((1-d,1+d), (-d,+d), **kwargs)
# ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)
#
# kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
# ax2.plot((-d,+d), (1-d,1+d), **kwargs)
# ax2.plot((-d,+d), (-d,+d), **kwargs)
# ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal




y_val_2 = 2

bax2 = brokenaxes(xlims=((-10000, -5500), (-500, 1000)), hspace=.05, despine=False)
for index, i in enumerate(epsilon_w182_5ma_dict['core']['epsilon_182w']):
    epsilon = epsilon_w182_5ma_dict['core']['epsilon_182w'][index]
    name = epsilon_w182_5ma_dict['core']['name'][index]
    if "IW" in name:
        bax2.scatter(epsilon, y_val_2, marker='D', s=80, color='blue')
    else:
        bax2.scatter(epsilon, y_val_2, marker='D', s=80, color='blue', facecolors="none")
    bax2.annotate(name, (epsilon + 100, y_val_2 - 0.4))
    y_val_2 += 2
for index, i in enumerate(epsilon_w182_5ma_dict['mantle']['epsilon_182w']):
    epsilon = epsilon_w182_5ma_dict['mantle']['epsilon_182w'][index]
    name = epsilon_w182_5ma_dict['mantle']['name'][index]
    if "IW" in name:
        bax2.scatter(epsilon, y_val_2, marker='s', s=80, color='green')
    else:
        bax2.scatter(epsilon, y_val_2, marker='s', s=80, color='green', facecolors="none")
    bax2.annotate(name, (epsilon + 100, y_val_2 - 0.4))
    y_val_2 += 2
for index, i in enumerate(epsilon_w182_5ma_dict['bulk']['epsilon_182w']):
    epsilon = epsilon_w182_5ma_dict['bulk']['epsilon_182w'][index]
    name = epsilon_w182_5ma_dict['bulk']['name'][index]
    if "IW" in name:
        bax2.scatter(epsilon, y_val_2, marker='o', s=80, color='red')
    else:
        bax2.scatter(epsilon, y_val_2, marker='o', s=80, color='red', facecolors="none")
    bax2.annotate(name, (epsilon + 100, y_val_2 - 0.4))
    y_val_2 += 2

bax2.grid(axis='x', which='major', ls='--', alpha=0.4)
axs = [i.yaxis.set_major_formatter(plt.NullFormatter()) for i in bax2.axs]
axs = [i.yaxis.set_ticks_position('none')  for i in bax2.axs]

bax2.set_title("$\epsilon$$^{182}$W at Time $t=5$ Ma")
bax2.set_xlabel("$\epsilon$$^{182}$W")
bax2.big_ax.xaxis.labelpad = 40






# 100 Ma below

# y_val_2 = 2
# 
# bax2 = brokenaxes(xlims=((-9900, -9800), (-1800, -1400), (-10, 115)), hspace=.0005, despine=False, d=0.01)
# for index, i in enumerate(epsilon_w182_100ma_dict['core']['epsilon_182w']):
#     epsilon = epsilon_w182_100ma_dict['core']['epsilon_182w'][index]
#     name = epsilon_w182_100ma_dict['core']['name'][index]
#     if "IW" in name:
#         bax2.scatter(epsilon, y_val_2, marker='D', s=80, color='blue')
#     else:
#         bax2.scatter(epsilon, y_val_2, marker='D', s=80, color='blue', facecolors="none")
#     bax2.annotate(name, (epsilon + 10, y_val_2 - 0.4))
#     y_val_2 += 2
# for index, i in enumerate(epsilon_w182_100ma_dict['mantle']['epsilon_182w']):
#     epsilon = epsilon_w182_100ma_dict['mantle']['epsilon_182w'][index]
#     name = epsilon_w182_100ma_dict['mantle']['name'][index]
#     if "IW" in name:
#         bax2.scatter(epsilon, y_val_2, marker='s', s=80, color='green')
#     else:
#         bax2.scatter(epsilon, y_val_2, marker='s', s=80, color='green', facecolors="none")
#     bax2.annotate(name, (epsilon + 10, y_val_2 - 0.4))
#     y_val_2 += 2
# for index, i in enumerate(epsilon_w182_100ma_dict['bulk']['epsilon_182w']):
#     epsilon = epsilon_w182_100ma_dict['bulk']['epsilon_182w'][index]
#     name = epsilon_w182_100ma_dict['bulk']['name'][index]
#     if "IW" in name:
#         bax2.scatter(epsilon, y_val_2, marker='o', s=80, color='red')
#     else:
#         bax2.scatter(epsilon, y_val_2, marker='o', s=80, color='red', facecolors="none")
#     bax2.annotate(name, (epsilon + 10, y_val_2 - 0.4))
#     y_val_2 += 2
# 
# bax2.grid(axis='x', which='major', ls='--', alpha=0.4)
# axs = [i.yaxis.set_major_formatter(plt.NullFormatter()) for i in bax2.axs]
# axs = [i.yaxis.set_ticks_position('none') for i in bax2.axs]
# 
# bax2.set_title("$\epsilon$$^{182}$W at Time $t=100$ Ma")
# bax2.set_xlabel("$\epsilon$$^{182}$W")
# bax2.big_ax.xaxis.labelpad = 40





# 5 Ma below

# y_val_2 = 2
#
# bax2 = brokenaxes(xlims=((-9900, -9500), (-7400, -6400), (-10, 200)), hspace=.0005, despine=False, d=0.01)
# for index, i in enumerate(epsilon_w182_5ma_dict['core']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['core']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['core']['name'][index]
#     bax2.scatter(epsilon, y_val_2, marker='^', s=80, color='blue')
#     bax2.annotate(name, (epsilon + 10, y_val_2 - 0.4))
#     y_val_2 += 2
# for index, i in enumerate(epsilon_w182_5ma_dict['mantle']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['mantle']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['mantle']['name'][index]
#     bax2.scatter(epsilon, y_val_2, marker='s', s=80, color='green')
#     bax2.annotate(name, (epsilon + 10, y_val_2 - 0.4))
#     y_val_2 += 2
# for index, i in enumerate(epsilon_w182_5ma_dict['bulk']['epsilon_182w']):
#     epsilon = epsilon_w182_5ma_dict['bulk']['epsilon_182w'][index]
#     name = epsilon_w182_5ma_dict['bulk']['name'][index]
#     bax2.scatter(epsilon, y_val_2, marker='o', s=80, color='red')
#     bax2.annotate(name, (epsilon + 10, y_val_2 - 0.4))
#     y_val_2 += 2
#
# bax2.grid(axis='x', which='major', ls='--', alpha=0.4)
# axs = [i.yaxis.set_major_formatter(plt.NullFormatter()) for i in bax2.axs]
# axs = [i.yaxis.set_ticks_position('none') for i in bax2.axs]
#
# bax2.set_title("$\epsilon$$^{182}$W at Time $t=5 Ma$")
# bax2.set_xlabel("$\epsilon$$^{182}$W")

plt.show()



