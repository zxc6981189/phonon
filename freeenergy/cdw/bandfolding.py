# This code aim to fold the band structure into the bigger supercell
import numpy as np
from lib import read_poscar

# the simple case
# cubic unit cell


def fold(length, kpt):
    """

    args:
        length(float): the cubic kspace length
        kpt(np.array): k point position in cartesian coordinate
    """
    return kpt-length*np.round(kpt/length)


def kspace(lattice):
    """
    transform real space lattice into k space

    args:
        lattice(np.array): lattice matrix in 3 by 3 array

    return:
        3 by 3 np.array in k space lattice
    """
    a1 = lattice[0]
    a2 = lattice[1]
    a3 = lattice[2]

    volume = np.dot(a1, np.cross(a2, a3))
    factor = 1/volume

    b1 = factor*np.cross(a2, a3)
    b2 = factor*np.cross(a3, a1)
    b3 = factor*np.cross(a1, a2)

    return np.array([b1, b2, b3])


def kpt(kmesh, klattice):
    """
    generate the k grid

    args:
        kmesh(np.array): k grid in 1 by 3 array
        lattice_constant(np.array): lattice constant 1 by 3 array

    return:
        kpoints list in direct coordinate
    """

    # remove the end point of liner space function
    nkpt = (kmesh[0]-1)*(kmesh[1]-1)*(kmesh[2]-1)

    # create the empty array to store the all kpoints
    kpoints = np.zeros((nkpt, 3))

    # length of lattice constant in k space
    b1 = np.linalg.norm(klattice[0])
    b2 = np.linalg.norm(klattice[1])
    b3 = np.linalg.norm(klattice[2])
    klattice_constant = np.array([b1, b2, b3])

    # linear space the point in x y z direction
    kx = np.linspace(-0.5, 0.5, kmesh[0])
    ky = np.linspace(-0.5, 0.5, kmesh[1])
    kz = np.linspace(-0.5, 0.5, kmesh[2])

    k_index = 0
    for i in range(kmesh[0]-1):
        for j in range(kmesh[1]-1):
            for k in range(kmesh[2]-1):
                kpoints[k_index] = [kx[i], ky[j], kz[k]]
                k_index += 1

    return kpoints*klattice_constant


if __name__ == '__main__':
    # lattice_constant = np.array([1, 1, 1])
    # lattice = np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    lattice = np.array(
        [[7.0747701492385300, 0, 0],
         [-3.5373850746192939, 6.126930675176626, 0],
            [0, 0, 12.8589043717472045]]
    )
    klattice = kspace(lattice)
    kpoints = kpt(np.array([25, 25, 13]), klattice)
