from math import pi, log, sqrt
import pandas as pd
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt
from decimal import Decimal

plt.rcParams.update({'font.size': 16})

def calcAge(half_life, ratio):
    decay_const = log(0.5) / half_life
    t = (1 / decay_const) * log(ratio)
    return t

ratios = list(np.arange(0.01, 1 + 0.01, 0.01))
half_life = 8.9 * (10**6)
ages = [calcAge(half_life=half_life, ratio=i) for i in ratios]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(ratios, [i / (10**6) for i in ages], linewidth=2.0, color='black')
ax1.set_title("Kleine et al. 2004 Core Age Sensitivity Modeling")
ax1.set_xlabel("$^{182}$Hf/$^{180}$Hf")
ax1.set_ylabel("Core Age (Ma)")
ax1.grid()


plt.show()