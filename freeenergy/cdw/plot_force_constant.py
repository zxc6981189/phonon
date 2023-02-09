import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lib import read_force_constant

fc_1 = read_force_constant(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-highk/13step/DP-0.01/FORCE_CONSTANTS')
fc_15 = read_force_constant(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-highk/13step/DP-0.15/FORCE_CONSTANTS')
fc_17 = read_force_constant(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-highk/13step/DP-0.17/FORCE_CONSTANTS')
fc_18 = read_force_constant(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-highk/13step/DP-0.18/FORCE_CONSTANTS')
fc_20 = read_force_constant(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-highk/13step/DP-0.20/FORCE_CONSTANTS')
fc_25 = read_force_constant(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-highk/13step/DP-0.25/FORCE_CONSTANTS')

# force constant tensor
# natom, natom, 3, 3
fc1 = fc_1[[0, 8, 16]]
fc15 = fc_15[[0, 8, 16]]
fc17 = fc_17[[0, 8, 16]]
fc18 = fc_18[[0, 8, 16]]
fc20 = fc_20[[0, 8, 16]]
fc25 = fc_25[[0, 8, 16]]

# percentage
# pfc1 = np.abs(fc1-fc1)/fc1
# pfc15 = np.abs(fc1-fc15)/fc1
# pfc17 = np.abs(fc1-fc17)/fc1
# pfc18 = np.abs(fc1-fc18)/fc1
# pfc20 = np.abs(fc1-fc20)/fc1
# pfc25 = np.abs(fc1-fc25)/fc1

# flatten the np array
rfc1 = fc1.reshape(-1)
rfc15 = fc15.reshape(-1)
rfc17 = fc17.reshape(-1)
rfc18 = fc18.reshape(-1)
rfc20 = fc20.reshape(-1)
rfc25 = fc25.reshape(-1)

x1 = np.repeat(0.01, 3*24)
# x15 = np.repeat(0.15, 576*9)
# x17 = np.repeat(0.17, 576*9)
# x18 = np.repeat(0.18, 576*9)
# x20 = np.repeat(0.20, 576*9)
# x25 = np.repeat(0.25, 576*9)

all_fc = [rfc1, rfc15, rfc17, rfc18, rfc20, rfc25]
all_fc = np.array(all_fc)

fig, ax = plt.subplots()

ax.plot(all_fc)
ax.set_xlabel('displacement data number')
ax.set_ylabel('force constant')

list_xaxis = [0.01, 0.15, 0.17, 0.18, 0.20, 0.25]
ax.set_xticks([0, 1, 2, 3, 4, 5])
ax.set_xticklabels(list_xaxis)
