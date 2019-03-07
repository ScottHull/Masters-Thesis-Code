import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, pi, sqrt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def iterReverseD(obj_concs, cell_concs, index=0, iterReverseDList=[]):

    if index < len(list(obj_concs)):
        obj = list(obj_concs)[index]
        cell_concs_range = list(cell_concs)[0:index + 1]
        avg_cell_concs_range = sum(cell_concs_range) / (len(cell_concs_range))
        avg_D = obj / avg_cell_concs_range
        iterReverseDList.append(avg_D)
        return iterReverseD(obj_concs=obj_concs, cell_concs=cell_concs, index=(index + 1), iterReverseDList=iterReverseDList)
    else:
        return iterReverseDList

def recalcConc(old_concentrations, old_volume, geom_type='box', radius=None, height=None, width=None, length=None):
    if geom_type == 'sphere' or geom_type == 's':
        volume = (4/3) * pi * (radius**3)
        new_conc = [(i * old_volume) / volume for i in list(old_concentrations)]
        return new_conc
    else:
        volume = length * width * height
        new_conc = [(i * old_volume) / volume for i in list(old_concentrations)]
        return new_conc

def calcDiffusionLength(chem_diffusivity, droplet_radius, settling_velocity):
    l = sqrt((2 * chem_diffusivity * droplet_radius) / settling_velocity)
    return l

def meltLengthWidth(diff_length, droplet_radius):
    length_width = (2 * droplet_radius) + (2 * diff_length)
    return length_width

df_fO2_08 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-0.8.csv")
df_fO2_110 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-1.10.csv")
df_fO2_225 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-2.25.csv")
df_fO2_245 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-2.45.csv")

D_fO2_08 = list(df_fO2_08['D'])[:-3]
depth_fO2_08 = list(df_fO2_08['z-depth'] / 1000)[:-3]
cell_temperatures_fO2_08 = df_fO2_08['cell_temperature']
cell_pressures_fO2_08 = df_fO2_08['cell_pressure']
obj_conc_fO2_08 = df_fO2_08['object_conc']
cell_conc_fO2_08 = df_fO2_08['cell_conc']
cell_moles_fO2_08 = list(df_fO2_08['cell_moles'])[:-3]
object_moles_fO2_08 = list(df_fO2_08['object_moles'])[:-3]
reverse_D_fO2_08 = iterReverseD(obj_concs=obj_conc_fO2_08, cell_concs=cell_conc_fO2_08, index=0, iterReverseDList=[])[:-3]

D_fO2_110 = list(df_fO2_110['D'])[:-3]
depth_fO2_110 = list(df_fO2_110['z-depth'] / 1000)[:-3]
cell_temperatures_fO2_110 = df_fO2_110['cell_temperature']
cell_pressures_fO2_110 = df_fO2_110['cell_pressure']
obj_conc_fO2_110 = df_fO2_110['object_conc']
cell_conc_fO2_110 = df_fO2_110['cell_conc']
cell_moles_fO2_110 = list(df_fO2_110['cell_moles'])[:-3]
object_moles_fO2_110 = list(df_fO2_110['object_moles'])[:-3]
reverse_D_fO2_110 = iterReverseD(obj_concs=obj_conc_fO2_110, cell_concs=cell_conc_fO2_110, index=0, iterReverseDList=[])[:-3]

D_fO2_225 = list(df_fO2_225['D'])[:-3]
depth_fO2_225 = list(df_fO2_225['z-depth'] / 1000)[:-3]
cell_temperatures_fO2_225 = df_fO2_225['cell_temperature']
cell_pressures_fO2_225 = df_fO2_225['cell_pressure']
obj_conc_fO2_225 = df_fO2_225['object_conc']
cell_conc_fO2_225 = df_fO2_225['cell_conc']
cell_moles_fO2_225 = list(df_fO2_225['cell_moles'])[:-3]
object_moles_fO2_225 = list(df_fO2_225['object_moles'])[:-3]
reverse_D_fO2_225 = iterReverseD(obj_concs=obj_conc_fO2_225, cell_concs=cell_conc_fO2_225, index=0, iterReverseDList=[])[:-3]

D_fO2_245 = list(df_fO2_245['D'])[:-3]
depth_fO2_245 = list(df_fO2_245['z-depth'] / 1000)[:-3]
cell_temperatures_fO2_245 = df_fO2_245['cell_temperature']
cell_pressures_fO2_245 = df_fO2_245['cell_pressure']
obj_conc_fO2_245 = df_fO2_245['object_conc']
cell_conc_fO2_245 = df_fO2_245['cell_conc']
cell_moles_fO2_245 = list(df_fO2_245['cell_moles'])[:-3]
object_moles_fO2_245 = list(df_fO2_245['object_moles'])[:-3]
reverse_D_fO2_245 = iterReverseD(obj_concs=obj_conc_fO2_245, cell_concs=cell_conc_fO2_245, index=0, iterReverseDList=[])[:-3]

diff_length = calcDiffusionLength(chem_diffusivity=10**-8, settling_velocity=0.25, droplet_radius=0.0185)
length = meltLengthWidth(droplet_radius=0.0185, diff_length=diff_length)
width = length

new_cell_conc_fO2_08 = recalcConc(old_concentrations=cell_conc_fO2_08, old_volume=400**3, length=length,
                                  width=width, height=400)
new_cell_conc_fO2_110 = recalcConc(old_concentrations=cell_conc_fO2_110, old_volume=400**3, length=length,
                                  width=width, height=400)
new_cell_conc_fO2_225 = recalcConc(old_concentrations=cell_conc_fO2_225, old_volume=400**3, length=length,
                                  width=width, height=400)
new_cell_conc_fO2_245 = recalcConc(old_concentrations=cell_conc_fO2_245, old_volume=400**3, length=length,
                                  width=width, height=400)


new_reverse_D_fO2_08 = iterReverseD(obj_concs=obj_conc_fO2_08, cell_concs=new_cell_conc_fO2_08, index=0, iterReverseDList=[])[:-3]
new_reverse_D_fO2_110 = iterReverseD(obj_concs=obj_conc_fO2_110, cell_concs=new_cell_conc_fO2_110, index=0, iterReverseDList=[])[:-3]
new_reverse_D_fO2_225 = iterReverseD(obj_concs=obj_conc_fO2_225, cell_concs=new_cell_conc_fO2_225, index=0, iterReverseDList=[])[:-3]
new_reverse_D_fO2_245 = iterReverseD(obj_concs=obj_conc_fO2_245, cell_concs=new_cell_conc_fO2_245, index=0, iterReverseDList=[])[:-3]

print(new_cell_conc_fO2_08)

fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax1.plot(depth_fO2_08, D_fO2_08, label='fO2 = IW-0.8', color='red')
ax1.plot(depth_fO2_08, new_reverse_D_fO2_08, label='Reversed fO2 = IW-0.8', color='red', linestyle="--")
ax1.plot(depth_fO2_110, D_fO2_110, label='fO2 = IW-1.10', color='blue')
ax1.plot(depth_fO2_110, new_reverse_D_fO2_110, label='Reversed fO2 = IW-1.10', color='blue', linestyle="--")
ax1_2 = fig1.add_subplot(212)
ax1_2.plot(depth_fO2_08, [i - j for i, j in zip(D_fO2_08, new_reverse_D_fO2_08)],
           label="Reversed - Distributed (fO2 = IW-0.8)", linestyle="--", color='red')
ax1_2.plot(depth_fO2_110, [i - j for i, j in zip(D_fO2_110, new_reverse_D_fO2_110)],
           label="Reversed - Distributed (fO2 = IW-1.10)", linestyle="--", color='blue')
ax1.grid()
ax1.fill_between(depth_fO2_08, D_fO2_08, D_fO2_110, color='red', alpha=0.2)
ax1_2.fill_between(depth_fO2_225, [i - j for i, j in zip(D_fO2_08, new_reverse_D_fO2_08)],
                   [i - j for i, j in zip(D_fO2_110, new_reverse_D_fO2_110)], color='red', alpha=0.2)
ax1.set_title("Iron Meteorite fO2 Model")
ax1.legend(loc='center right')
ax1_2.legend(loc='center right')
ax1_2.set_xlabel("Depth (km)")
ax1.set_ylabel("Bulk D ('Reversed')")
ax1_2.set_ylabel("Reversed - Distributed")
ax1_2.grid()


fig2 = plt.figure()
ax2 = fig2.add_subplot(211)
ax2.plot(depth_fO2_225, D_fO2_225, label='fO2 = IW-2.25', color='red')
ax2.plot(depth_fO2_225, new_reverse_D_fO2_225, label='Reversed fO2 = IW-2.25', color='red', linestyle="--")
ax2.plot(depth_fO2_245, D_fO2_245, label='fO2 = IW-2.45', color='blue')
ax2.plot(depth_fO2_245, new_reverse_D_fO2_245, label='Reversed fO2 = IW-2.45', color='blue', linestyle="--")
ax2_2 = fig2.add_subplot(212)
ax2_2.plot(depth_fO2_225, [i - j for i, j in zip(D_fO2_225, new_reverse_D_fO2_225)],
           label="Reversed - Distributed (fO2 = IW-2.25)", linestyle="--", color='red')
ax2_2.plot(depth_fO2_245, [i - j for i, j in zip(D_fO2_245, new_reverse_D_fO2_245)],
           label="Reversed - Distributed (fO2 = IW-2.45)", linestyle="--", color='blue')
ax2.grid()
ax2.fill_between(depth_fO2_225, D_fO2_225, D_fO2_245, color='red', alpha=0.2)
ax2_2.fill_between(depth_fO2_225, [i - j for i, j in zip(D_fO2_225, new_reverse_D_fO2_225)],
                   [i - j for i, j in zip(D_fO2_245, new_reverse_D_fO2_245)], color='red', alpha=0.2)
ax2.legend(loc='center right')
ax2_2.legend(loc='center right')
ax2.set_title("Vesta fO2 Model")
ax2_2.set_xlabel("Depth (km)")
ax2.set_ylabel("Bulk D ('Reversed')")
ax2_2.set_ylabel("Reversed - Distributed")
ax2_2.grid()



plt.show()