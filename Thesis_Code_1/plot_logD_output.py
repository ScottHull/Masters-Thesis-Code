import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd

def calcBulkConc(conc_cells, conc_objs):
    bulk_conc = [obj_conc + cell_conc for obj_conc, cell_conc in zip(conc_objs, conc_cells)]
    return bulk_conc

def calcAdjBulkConc(massfrac_silicate, massfrac_metal, conc_cells, conc_objs):
    adj_bulk_conc = [((massfrac_metal * obj_conc) + (massfrac_silicate * cell_conc))
                     for obj_conc, cell_conc in zip(conc_objs, conc_cells)]
    return adj_bulk_conc


def calcAdjD(massfrac_silicate, massfrac_metal, obj_concs, adj_bulk_concs):
    adj_D = [((massfrac_silicate * obj_conc) / adj_bulk_conc) - (massfrac_silicate / massfrac_metal) for
             obj_conc, adj_bulk_conc in zip(obj_concs, adj_bulk_concs)]
    return adj_D

def reverseD(obj_concs, cell_concs):
    final_obj_conc = list(obj_concs)[-1]
    # print("~~~  {} ~~~".format(final_obj_conc))
    rev_cell_concs = [final_obj_conc / i for i in cell_concs]
    return rev_cell_concs

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







df_fO2_245 = pd.read_csv("single_droplet/D_coeffs_fO2=-2.45.csv")
df_fO2_225 = pd.read_csv("single_droplet/D_coeffs_fO2=-2.25.csv")
df_fO2_110 = pd.read_csv("single_droplet/D_coeffs_fO2=-1.10.csv")
df_fO2_08 = pd.read_csv("single_droplet/D_coeffs_fO2=-0.8.csv")

massfrac_silicate = 0.68
massfrac_metal = 0.32

D_fO2_245 = df_fO2_245['D']
depth_fO2_245 = df_fO2_245['z-depth'] / 1000
cell_temperatures_fO2_245 = df_fO2_245['cell_temperature']
cell_pressures_fO2_245 = df_fO2_245['cell_pressure']
obj_conc_fO2_245 = df_fO2_245['object_conc']
cell_conc_fO2_245 = df_fO2_245['cell_conc']
bulk_conc_fO2_245 = calcBulkConc(conc_cells=cell_conc_fO2_245, conc_objs=obj_conc_fO2_245)
adj_bulk_conc_fO2_245 = calcAdjBulkConc(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
                                        conc_cells=cell_conc_fO2_245, conc_objs=obj_conc_fO2_245)
adj_D_fO2_245 = calcAdjD(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
                                        adj_bulk_concs=adj_bulk_conc_fO2_245, obj_concs=obj_conc_fO2_245)
rev_cell_conc_fO2_245 = reverseD(obj_concs=obj_conc_fO2_245, cell_concs=cell_conc_fO2_245)
iter_rev_cell_conc_fO2_245 = iterReverseD(obj_concs=obj_conc_fO2_245, cell_concs=cell_conc_fO2_245)

D_fO2_225 = df_fO2_225['D']
depth_fO2_225 = df_fO2_225['z-depth'] / 1000
cell_temperatures_fO2_225 = df_fO2_225['cell_temperature']
cell_pressures_fO2_225 = df_fO2_225['cell_pressure']
obj_conc_fO2_225 = df_fO2_225['object_conc']
cell_conc_fO2_225 = df_fO2_225['cell_conc']
bulk_conc_fO2_225 = calcBulkConc(conc_cells=cell_conc_fO2_225, conc_objs=obj_conc_fO2_225)
adj_bulk_conc_fO2_225 = calcAdjBulkConc(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
                                        conc_cells=cell_conc_fO2_225, conc_objs=obj_conc_fO2_225)
adj_D_fO2_225 = calcAdjD(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
                                        adj_bulk_concs=adj_bulk_conc_fO2_225, obj_concs=obj_conc_fO2_225)
rev_cell_conc_fO2_225 = reverseD(obj_concs=obj_conc_fO2_225, cell_concs=cell_conc_fO2_225)
# iter_rev_cell_conc_fO2_225 = iterReverseD(obj_concs=obj_conc_fO2_225, cell_concs=cell_conc_fO2_225)

#
# D_fO2_110 = df_fO2_110['D']
# depth_fO2_110 = df_fO2_110['z-depth'] / 1000
# cell_temperatures_fO2_110 = df_fO2_110['cell_temperature']
# cell_pressures_fO2_110 = df_fO2_110['cell_pressure']
# obj_conc_fO2_110 = df_fO2_110['object_conc']
# cell_conc_fO2_110 = df_fO2_110['cell_conc']
# bulk_conc_fO2_110 = calcBulkConc(conc_cells=cell_conc_fO2_110, conc_objs=obj_conc_fO2_110)
# adj_bulk_conc_fO2_110 = calcAdjBulkConc(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
#                                         conc_cells=cell_conc_fO2_110, conc_objs=obj_conc_fO2_110)
# adj_D_fO2_110 = calcAdjD(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
#                                         adj_bulk_concs=adj_bulk_conc_fO2_110, obj_concs=obj_conc_fO2_110)
# rev_cell_conc_fO2_110 = reverseD(obj_concs=obj_conc_fO2_110, cell_concs=cell_conc_fO2_110)
#
# D_fO2_08 = df_fO2_08['D']
# depth_fO2_08 = df_fO2_08['z-depth'] / 1000
# cell_temperatures_fO2_08 = df_fO2_08['cell_temperature']
# cell_pressures_fO2_08 = df_fO2_08['cell_pressure']
# obj_conc_fO2_08 = df_fO2_08['object_conc']
# cell_conc_fO2_08 = df_fO2_08['cell_conc']
# bulk_conc_fO2_08 = calcBulkConc(conc_cells=cell_conc_fO2_08, conc_objs=obj_conc_fO2_08)
# adj_bulk_conc_fO2_08 = calcAdjBulkConc(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
#                                         conc_cells=cell_conc_fO2_08, conc_objs=obj_conc_fO2_08)
# adj_D_fO2_08 = calcAdjD(massfrac_silicate=massfrac_silicate, massfrac_metal=massfrac_metal,
#                                         adj_bulk_concs=adj_bulk_conc_fO2_08, obj_concs=obj_conc_fO2_08)
# rev_cell_conc_fO2_08 = reverseD(obj_concs=obj_conc_fO2_08, cell_concs=cell_conc_fO2_08)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
# ax1.plot(depth_fO2_245, D_fO2_245, linewidth=2.0, label="fO2 IW -2.45", linestyle="--")
# ax1.plot(depth_fO2_225, D_fO2_225, linewidth=2.0, label="fO2 IW -2.25", linestyle="--")
ax1.plot(depth_fO2_245, [i - j for i,j in zip(iter_rev_cell_conc_fO2_245, D_fO2_245)], linewidth=2.0, label='Iter Rev D IW -2.45')
# ax1.plot(depth_fO2_245, rev_cell_conc_fO2_245, linewidth=2.0, label="adj fO2 IW -2.45")
# ax1.fill_between(depth_fO2_225, D_fO2_225, D_fO2_245, color='red', alpha=0.2)
ax1.set_xlabel("Depth (km)")
ax1.set_ylabel("Partition Coefficient, D")
ax1.set_title("Vestian W Partition Coefficients (Steenstra et al. 2016 fO2 Range)")
ax1.grid()
ax1.legend(loc='center right')

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# ax2.plot(depth_fO2_110, D_fO2_110, linewidth=2.0, label="fO2 IW -1.10", linestyle="--")
# ax2.plot(depth_fO2_08, D_fO2_08, linewidth=2.0, label="fO2 IW -0.80", linestyle="--")
# ax2.fill_between(depth_fO2_110, D_fO2_08, D_fO2_110, color='red', alpha=0.2)
# ax2.set_xlabel("Depth (km)")
# ax2.set_ylabel("Partition Coefficient, D")
# ax2.set_title("Vestian W Partition Coefficients (Bonnand et al. 2018 fO2 Range)")
# ax2.grid()
# ax2.legend(loc='upper left')

plt.show()