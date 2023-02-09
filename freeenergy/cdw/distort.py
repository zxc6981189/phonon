#! /Users/liangying/Documents/Code/FreeEnergy/.venv/bin/python
from lib import write_poscar
from lib import read_poscar
from lib import vector_displacement


def distorted_poscar(nondistorted, distorted, percentage):
    """
    args:
        nondistorted(str): filename of nondistorted poscar
        distorted(str): filename of distorted poscar
        percentage(float): the distorted percentage
        e.g. 0.1, shift atoms from nondistorted to distorted
        in 10% displacement
    return:
        (np.array): distorted cdw crystal structure
    """

    # nondistorted poscar
    ndispos = read_poscar(nondistorted, True)

    # position changed in between two poscar
    displacement = vector_displacement(nondistorted, distorted)

    # move atoms in certain percentage of displacement
    # in cartesian coordinate
    percentage_position = ndispos['position']+displacement*percentage

    # replace position in poscar dict
    ndispos['position'] = percentage_position

    # output the percentage distorted poscar
    # write_poscar('POSCAR-dis', ndispos, True)
    return percentage_position
