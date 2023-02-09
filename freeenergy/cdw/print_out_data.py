#! /Users/liangying/Documents/Code/FreeEnergy/.venv/bin/python
from lib import hdf5_parse


def print_frequency(data, segment, point, nband):
    """
    This code is to print out the frequency at specific point

    args:
        data(array): a phonon frequency numpy array object
                     output from hdf5_parse
        segment(int): the band interval start from 0
        point(int): the point index in the interval start from 0
        nband(int): number of the band start from 0
    """

    print(data[segment, point, nband])


data = hdf5_parse('band.hdf5')

# Gamma point
print('Gamma')
print_frequency(data['frequency'], 3, 0, 0)
# A point
print('A')
print_frequency(data['frequency'], 3, -1, 0)
# Gamma point
print('L')
print_frequency(data['frequency'], 4, -1, 0)
