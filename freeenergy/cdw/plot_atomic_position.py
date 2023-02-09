import numpy as np
from lib import read_poscar
from lib import supercell
import matplotlib.pyplot as plt

# read poscar

distort_poscar = read_poscar(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-2x2x2/ISIF2-GGA-sym-highk/5step/1-Sc/POSCAR')
normal_poscar = read_poscar(
    '/Users/liangying/Documents/Server/T3/Work/TiSe2/bulk-1x1x1/13step/1-Sc/POSCAR')

# read the data out
distort_lattice = distort_poscar['lattice']
normal_lattice = normal_poscar['lattice']

distort_position = distort_poscar['position']
normal_position = normal_poscar['position']

distort_supercell_lattice, distort_supercell_position = supercell(
    [1, 1, 1], distort_lattice, distort_position)
normal_supercell_lattice, normal_supercell_position = supercell(
    [1, 1, 1], normal_lattice, normal_position)

cart_distort_position = np.dot(
    distort_supercell_position, distort_supercell_lattice)
cart_normal_position = np.dot(
    normal_supercell_position, normal_supercell_lattice)


plt.scatter(cart_distort_position[:, 0],
            cart_distort_position[:, 1],
            c='b')
plt.scatter(cart_normal_position[:, 0],
            cart_normal_position[:, 1],
            c='r')

plt.show()
