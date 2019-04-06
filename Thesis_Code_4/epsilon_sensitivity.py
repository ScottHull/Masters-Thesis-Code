import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, sqrt, pi, log10
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rcParams.update({'font.size': 16})

def epsilon182W(sample_182W_184W, standard_182W_184W):
    epsilon = ((sample_182W_184W / standard_182W_184W) - 1) * (10**4)
    return epsilon


terrestrial_standard_182W_184W = 0.864680
sample_182W_184W = np.arange(terrestrial_standard_182W_184W - 0.01, terrestrial_standard_182W_184W + 0.01, 0.001)
print(sample_182W_184W)

epsilons = [epsilon182W(sample_182W_184W=i, standard_182W_184W=terrestrial_standard_182W_184W) for i in sample_182W_184W]

print(epsilons)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(sample_182W_184W, epsilons, linewidth=2.0, color='black')
ax1.axvline(terrestrial_standard_182W_184W, linewidth=2.0, linestyle="--", color='black',
            label='Terrestrial Standard (Kleine et al. 2004)')
ax1.set_xlabel("$^{182}$W/$^{184}$W")
ax1.set_ylabel("$\epsilon^{182}$W")
ax1.set_title("Sensitivity of $\epsilon^{182}$W to Error")
ax1.grid()
ax1.legend(loc='upper left')

plt.show()
