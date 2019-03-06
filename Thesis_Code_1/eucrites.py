import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, pi
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

def calcObjMoleFlux(object_moles, index=0, object_mole_flux=[]):
    if index == 0:
        object_mole_flux.append(0)
        return calcObjMoleFlux(object_moles=object_moles, index=(index + 1), object_mole_flux=object_mole_flux)
    elif index < len(object_moles):
        object_mole_flux.append(object_moles[index] - object_moles[index - 1])
        return calcObjMoleFlux(object_moles=object_moles, index=(index + 1), object_mole_flux=object_mole_flux)
    else:
        return object_mole_flux

def calcNumTotalDroplets(core_radius, droplet_radius):
    core_volume = (4/3) * pi * (core_radius**3)
    droplet_volume = (4/3) * pi * (droplet_radius**3)
    num_droplets = core_volume / droplet_volume
    return num_droplets


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


fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax1.plot(depth_fO2_08, D_fO2_08, label='fO2 = IW-0.8', color='red')
ax1.plot(depth_fO2_08, reverse_D_fO2_08, label='Reversed fO2 = IW-0.8', color='red', linestyle="--")
ax1.plot(depth_fO2_110, D_fO2_110, label='fO2 = IW-1.10', color='blue')
ax1.plot(depth_fO2_110, reverse_D_fO2_110, label='Reversed fO2 = IW-1.10', color='blue', linestyle="--")
ax1_2 = fig1.add_subplot(212)
ax1_2.plot(depth_fO2_08, [i - j for i, j in zip(D_fO2_08, reverse_D_fO2_08)],
           label="Reversed - Distributed (fO2 = IW-0.8)", linestyle="--", color='red')
ax1_2.plot(depth_fO2_110, [i - j for i, j in zip(D_fO2_110, reverse_D_fO2_110)],
           label="Reversed - Distributed (fO2 = IW-1.10)", linestyle="--", color='blue')
ax1.grid()
ax1.fill_between(depth_fO2_08, D_fO2_08, D_fO2_110, color='red', alpha=0.2)
ax1_2.fill_between(depth_fO2_225, [i - j for i, j in zip(D_fO2_08, reverse_D_fO2_08)],
                   [i - j for i, j in zip(D_fO2_110, reverse_D_fO2_110)], color='red', alpha=0.2)
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
ax2.plot(depth_fO2_225, reverse_D_fO2_225, label='Reversed fO2 = IW-2.25', color='red', linestyle="--")
ax2.plot(depth_fO2_245, D_fO2_245, label='fO2 = IW-2.45', color='blue')
ax2.plot(depth_fO2_245, reverse_D_fO2_245, label='Reversed fO2 = IW-2.45', color='blue', linestyle="--")
ax2_2 = fig2.add_subplot(212)
ax2_2.plot(depth_fO2_225, [i - j for i, j in zip(D_fO2_225, reverse_D_fO2_225)],
           label="Reversed - Distributed (fO2 = IW-2.25)", linestyle="--", color='red')
ax2_2.plot(depth_fO2_245, [i - j for i, j in zip(D_fO2_245, reverse_D_fO2_245)],
           label="Reversed - Distributed (fO2 = IW-2.45)", linestyle="--", color='blue')
ax2.grid()
ax2.fill_between(depth_fO2_225, D_fO2_225, D_fO2_245, color='red', alpha=0.2)
ax2_2.fill_between(depth_fO2_225, [i - j for i, j in zip(D_fO2_225, reverse_D_fO2_225)],
                   [i - j for i, j in zip(D_fO2_245, reverse_D_fO2_245)], color='red', alpha=0.2)
ax2.legend(loc='center right')
ax2_2.legend(loc='center right')
ax2.set_title("Vesta fO2 Model")
ax2_2.set_xlabel("Depth (km)")
ax2.set_ylabel("Bulk D ('Reversed')")
ax2_2.set_ylabel("Reversed - Distributed")
ax2_2.grid()

obj_mole_flux_08 = calcObjMoleFlux(object_moles=object_moles_fO2_08, index=0, object_mole_flux=[])
obj_mole_flux_110 = calcObjMoleFlux(object_moles=object_moles_fO2_110, index=0, object_mole_flux=[])
obj_mole_flux_225 = calcObjMoleFlux(object_moles=object_moles_fO2_225, index=0, object_mole_flux=[])
obj_mole_flux_245 = calcObjMoleFlux(object_moles=object_moles_fO2_245, index=0, object_mole_flux=[])

num_dropelts_vesta = calcNumTotalDroplets(core_radius=113 * 1000, droplet_radius=0.0185)

adj_moles_fO2_08 = [i * num_dropelts_vesta for i in object_moles_fO2_08]
adj_moles_fO2_110 = [i * num_dropelts_vesta for i in object_moles_fO2_110]
adj_moles_fO2_225 = [i * num_dropelts_vesta for i in object_moles_fO2_225]
adj_moles_fO2_245 = [i * num_dropelts_vesta for i in object_moles_fO2_245]
adj_obj_mole_flux_08 = calcObjMoleFlux(object_moles=adj_moles_fO2_08, index=0, object_mole_flux=[])
adj_obj_mole_flux_110 = calcObjMoleFlux(object_moles=adj_moles_fO2_110, index=0, object_mole_flux=[])
adj_obj_mole_flux_225 = calcObjMoleFlux(object_moles=adj_moles_fO2_225, index=0, object_mole_flux=[])
adj_obj_mole_flux_245 = calcObjMoleFlux(object_moles=adj_moles_fO2_245, index=0, object_mole_flux=[])


fig3 = plt.figure()
ax3 = fig3.add_subplot(211)
ax3_2 = fig3.add_subplot(212)
ax3.plot(depth_fO2_08, object_moles_fO2_08, color='red', label='fO2 = IW-0.8')
ax3.plot(depth_fO2_110, object_moles_fO2_110, color='blue', label='fO2 = IW-1.10')
ax3.plot(depth_fO2_225, object_moles_fO2_225, color='green', label='fO2 = IW-2.25')
ax3.plot(depth_fO2_245, object_moles_fO2_245, color='orange', label='fO2 = IW-2.45')
ax3_2.plot(depth_fO2_08, obj_mole_flux_08, color='red', label='fO2 = IW-0.8')
ax3_2.plot(depth_fO2_110, obj_mole_flux_110, color='blue', label='fO2 = IW-1.10')
ax3_2.plot(depth_fO2_225, obj_mole_flux_225, color='green', label='fO2 = IW-2.25')
ax3_2.plot(depth_fO2_245, obj_mole_flux_245, color='orange', label='fO2 = IW-2.45')
ax3.set_title("Single Droplet in the Vestian Magma Ocean")
ax3.grid()
ax3_2.grid()
ax3_2.set_xlabel("Depth (km)")
ax3.set_ylabel("Moles in Metal")
ax3_2.set_ylabel("Mole Flux into Metal")
ax3.legend(loc='center right')
ax3_2.legend(loc='center right')

fig4 = plt.figure()
ax4 = fig4.add_subplot(211)
ax4_2 = fig4.add_subplot(212)
ax4.plot(depth_fO2_08, adj_moles_fO2_08, color='red', label='fO2 = IW-0.8')
ax4.plot(depth_fO2_110, adj_moles_fO2_110, color='blue', label='fO2 = IW-1.10')
ax4.plot(depth_fO2_225, adj_moles_fO2_225, color='green', label='fO2 = IW-2.25')
ax4.plot(depth_fO2_245, adj_moles_fO2_245, color='orange', label='fO2 = IW-2.45')
ax4_2.plot(depth_fO2_08, adj_obj_mole_flux_08, color='red', label='fO2 = IW-0.8')
ax4_2.plot(depth_fO2_110, adj_obj_mole_flux_110, color='blue', label='fO2 = IW-1.10')
ax4_2.plot(depth_fO2_225, adj_obj_mole_flux_225, color='green', label='fO2 = IW-2.25')
ax4_2.plot(depth_fO2_245, adj_obj_mole_flux_245, color='orange', label='fO2 = IW-2.45')
ax4.set_title("Extrapolation to # of Droplets Reflective of Vesta's Core: {}".format(round(num_dropelts_vesta, 4)))
ax4.grid()
ax4_2.grid()
ax4_2.set_xlabel("Depth (km)")
ax4.set_ylabel("Moles in Metal")
ax4_2.set_ylabel("Mole Flux into Metal")
ax4.legend(loc='center right')
ax4_2.legend(loc='center right')


plt.show()