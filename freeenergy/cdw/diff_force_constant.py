#! /Users/liangying/Documents/Code/FreeEnergy/.venv/bin/python
# this program is to compare two force constant
from lib import read_force_constant
import numpy as np
import sys


def diff_force_constant(file_previous, file_current):
    """
    calculate the difference of force constant

    args:
        file_previous: previous self-consistent loop FORCE_CONTANTS file
        file_current: current self-consistent loop FORCE_CONTANTS file

    return:
        maxium and minium difference in the full force constant array
    """
    f1 = read_force_constant(file_previous)
    f2 = read_force_constant(file_current)
    diff = np.abs(f1-f2)
    return (np.amax(diff), np.mean(diff), np.nanmax(diff/f1))


if __name__ == '__main__':
    fc_max, fc_mean, fc_r_max = diff_force_constant(sys.argv[1], sys.argv[2])
    print('force constants absolute maxium: ', fc_max)
    print('force constants absolute mean: ', fc_mean)
    print('force constants relative maxium: ', fc_r_max)
