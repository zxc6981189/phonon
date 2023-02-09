#! /Users/liangying/Documents/Code/FreeEnergy/.venv/bin/python
import numpy as np
import time


def read_line(line: str, data_type='float'):
    """
    sub function to split the words into a list in float number

    args:
        line(str): line of string read from file

    return:
        list: the transformed float or int number list
    """

    if data_type == 'float':
        tmp_line = [float(x) for x in line.split()]
    elif data_type == 'int':
        tmp_line = [int(x) for x in line.split()]
    elif data_type == 'str':
        tmp_line = line.split()

    return tmp_line


def read_force_constant(filename: str):
    """
    read force constant from file

    args:
        filename(str): force constant file name (FORCE_CONSTANT)

    return:
        force constant in np.array (natom, natom, 3, 3)
    """

    with open(filename, 'r') as file:
        # line 1
        # number of atoms, number of atoms
        natom_a, natom_b = read_line(file.readline(), data_type='int')

        force_constant = np.zeros((natom_a, natom_b, 3, 3))

        for tmp in range(int(natom_a*natom_b)):
            index_natom_a, index_natom_b = read_line(
                file.readline(), data_type='int')

            # force constant tensor
            xx, xy, xz = read_line(file.readline())
            yx, yy, yz = read_line(file.readline())
            zx, zy, zz = read_line(file.readline())

            # the atomic index in the FORCE_CONSTANTS file is start from 1
            # Warning: the partial force constant file
            # cannot use the number of atom as an index
            force_constant[index_natom_a-1,
                           index_natom_b-1, 0, :] = [xx, xy, xz]
            force_constant[index_natom_a-1,
                           index_natom_b-1, 1, :] = [yx, yy, yz]
            force_constant[index_natom_a-1,
                           index_natom_b-1, 2, :] = [zx, zy, zz]
        return force_constant


def write_force_constant(FC, filename='FORCE_CONSTANTS_NEW'):
    """
    output the force constant from np array (FC)
    into a file named in filename variable

    args:
        FC(np.array): force constant numpy array
        filename(str): output file name (default: FORCE_CONSTANTS_NEW)

    return:
        output a force constant file
    """
    natom, _, _, _ = FC.shape

    with open(filename, 'w') as file:
        file.writelines('{:>4} {:>4}\n'.format(natom, natom))
        for col in range(natom):
            for row in range(natom):
                file.writelines('{} {}\n'.format(col+1, row+1))
                text = '    {fc[0]:18.15f}    {fc[1]:18.15f}    {fc[2]:18.15f}\n'
                file.writelines(text.format(fc=FC[col, row, 0]))
                file.writelines(text.format(fc=FC[col, row, 1]))
                file.writelines(text.format(fc=FC[col, row, 2]))


def write_force_constant_mix(previous_fc, current_fc, mixing=0.1):
    """
    output the force constant with mixing parameter

    args:
        previous_fc(str): previous force constant file
        current_fc(str): current force constant file

    return:
        output a mixing force constant
    """
    prev_fc = read_force_constant(previous_fc)
    curr_fc = read_force_constant(current_fc)

    mix_fc = (1-mixing)*prev_fc + mixing*curr_fc
    write_force_constant(mix_fc, filename='FORCE_CONSTANTS_MIX')


def direct_to_cart(lattice, position):
    """
    This function is to transformed
    the direct coordinate into cartesian coordinate

    args:
        lattice(np.array): the lattice matrix
        position(np.arrry): atomic positions in direct coordinate

    return:
        np.array: atomic positions in cartesian coordinate
    """
    cart = np.dot(position, lattice)

    return cart


def cart_to_direct(lattice, position):
    """
    This function is to transformed
    the cartesian coordinate into direct coordinate

    args:
        lattice(np.array): the lattice matrix
        position(np.arrry): atomic positions in direct coordinate

    return:
        np.array: atomic positions in cartesian coordinate
    """
    direct = np.dot(position, np.linalg.inv(lattice))

    return direct


def read_poscar(filename, cart=False):
    """
    read from the poscar
    args:
        filename(str): path of POSCAR
        cart(bol): transform atomic position into cartesian coordinate

    return:
        dict: dictionary include all the information of POSCAR
    """

    with open(filename, 'r') as file:
        # skip the name of system
        file.readline()

        # ratio
        ratio = read_line(file.readline())

        # lattice matrix
        lattice = np.zeros((3, 3), dtype='float64')
        for i in range(3):
            lattice[i] = read_line(file.readline())

        # atom type
        atom_type = read_line(file.readline(), data_type='str')

        # number of atom
        natom = read_line(file.readline(), data_type='int')

        # coordinate: direct or cart
        coordinate = read_line(file.readline(), data_type='str')

        total_natom = np.sum(natom)
        position = np.zeros((total_natom, 3))

        for i in range(total_natom):
            position[i] = read_line(file.readline())

        if cart:
            position = direct_to_cart(lattice, position)

    return {'ratio': ratio, 'lattice': lattice,
            'type': atom_type, 'natom': natom,
            'coordinate': coordinate, 'position': position}


def write_poscar(filename, poscar, cart=False):
    """
    write the poscar
    args:
        filename(str): output written poscar filename
        poscar(dict): poscar dictionary object
        cart(bol): default=False, if the position from
        poscar object is cartesian coordinate, use True

    return:
        output poscar file from certain poscar object
    """

    with open(filename, 'w') as file:
        # output time
        now = time.strftime("%a, %d %b %Y %H:%M:%S CST", time.localtime())

        # first line of file
        file.writelines('this poscar is output at '+now)
        file.writelines('\n')

        # ratio
        file.writelines('{:>19.14f}'.format(poscar['ratio'][0]))
        file.writelines('\n')

        # lattice matrix
        for i in range(3):
            file.writelines(' {lat[0]:22.16f}{lat[1]:22.16f}{lat[2]:22.16f}'.format(
                lat=poscar['lattice'][i, :].tolist()))
            file.writelines('\n')

        # atom type
        ntype = len(poscar['type'])
        for i in range(ntype):
            file.writelines('{:>5s}'.format(poscar['type'][i]))
        file.writelines('\n')

        # number of atom
        for i in range(ntype):
            file.writelines(' {:5d}'.format(poscar['natom'][i]))
        file.writelines('\n')

        # coordinate: direct or cart
        file.writelines('{:s}'.format(poscar['coordinate'][0]))
        file.writelines('\n')

        if cart:
            poscar['position'] = cart_to_direct(
                poscar['lattice'], poscar['position'])
        natom = np.sum(poscar['natom'])
        for i in range(natom):
            file.writelines(' {pos[0]:20.16f}{pos[1]:20.16f}{pos[2]:20.16f}'.format(
                pos=poscar['position'][i, :]))
            file.writelines('\n')


def supercell(size, lattice, position):
    """
    extend the unit cell to n by n by n super cell

    args:
        size(list): (list in length 3) desire extended size of supercell.
            e.g. [3, 3, 3]. The function will return all the atomic position
            in the supercell
        lattice(np.array): (3, 3) array. Unit cell lattice array
        position(np.array): atomic position in cartesian coordinate

    return:
        np.array: (n, 3) array. extended super cell atomic position
    """
    # number of atom
    natom = position.shape[0]

    # transform the coordinate

    position = direct_to_cart(lattice, position)

    # the factor to extend the unit cell
    ri, rj, rk = size

    # supercell lattice
    supercell_lattice = np.array(size)*lattice

    list_extend = []
    for i in range(ri):
        for j in range(rj):
            for k in range(rk):
                list_extend.append([i, j, k])

    # create a new list to save all the supercell atomic position
    supercell_position = []
    for factor in list_extend:

        # shifting vector
        shift = np.dot(factor, lattice)

        # shift the atom to the new position
        new_position = position+shift

        for index_atom in range(natom):
            supercell_position.append(new_position[index_atom].tolist())

    return supercell_lattice, supercell_position


def vector_displacement(poscar_normal, poscar_distort):
    """
    This function will return the displacement between two poscar position

    arg:
        poscar_normal(str): nondistorted poscar filename
        poscar_distort(srt): distorted poscar filename
    """

    position_normal = read_poscar(poscar_normal)['position']
    position_distort = read_poscar(poscar_distort)['position']

    lattice_normal = read_poscar(poscar_normal)['lattice']
    lattice_distort = read_poscar(poscar_distort)['lattice']

    position_normal = direct_to_cart(lattice_normal, position_normal)
    position_distort = direct_to_cart(lattice_distort, position_distort)

    displacement = position_distort-position_normal

    return displacement


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
    ndispos['position'] += displacement*percentage

    # output the percentage distorted poscar
    write_poscar('POSCAR-dis', ndispos, cart=True)


def hdf5_parse(filename='band.hdf5'):
    import h5py
    """
    parse the hdf5 file

    return:
        data: hdf5 data
    """
    data = h5py.File(filename)
    return data
