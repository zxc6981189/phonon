import numpy as np
from lib import hdf5_parse
from lib import read_poscar
from lib import cart_to_direct


def eigenvector_direct_coord(file_poscar, file_hdf5, cart=False):
    """
    this function is to convert the eigenvector list into the direct coordinate

    args:
        file_poscar(str): poscar file name
        file_hdf5(str): hdf5 data file name output from phonopy
    return:
        np.array: the eigenvector array (ninterval, npoint, nband, natom, 3)
        in direct coordinate
    """
    # read the poscar and hdf5 data
    poscar = read_poscar(file_poscar)
    data = hdf5_parse(file_hdf5)['eigenvector']

    # transform the array in to np.array
    data = np.array(data)

    # np.array: interval, npoint, natom*3dim, nband for each atoms
    # for example: data[0,0,0,0] is the first q point in the first interval
    # the first band of frist atom along x direction
    ninterval, npoint, neig, nband = data.shape

    # number of atom
    natom = int(neig/3)

    # rearrange the data
    for a in range(ninterval):
        for b in range(npoint):
            data[a, b] = np.transpose(data[a, b])

    # reshape the array
    data = data.reshape(ninterval, npoint, nband, natom, 3)

    # transform the coordinate
    if not cart:
        for a in range(ninterval):
            for b in range(npoint):
                for c in range(nband):
                    for d in range(natom):
                        tmp = cart_to_direct(
                            poscar['lattice'], data[a, b, c, d])
                        data[a, b, c, d, :] = tmp
    return data


def replace(poscar_vesta, data, interval, point, band, natom, scale=5):
    """
    This function will replace the content under VECTR and VECTT

    args:
        poscar_vesta(str): The POSCAR.vesta file name
        data(eigenvector): eigenvector object
                           output from eigenvector_direct_coord
        interval(int): the chose interval index
        point(int): the chose point index within the interval
        band(int): the chose band index
        natom(int): total number of atoms in the system
        scale(float): scale number. To enlarge the vector length

    return:
        output the replaced POSCAR.vesta file name in POSCAR_edit.vesta
    """

    # read all lines in the files
    with open(poscar_vesta, 'r') as file:
        # content is a list contant the lines one by one
        content = file.readlines()

    # find the key word 'VECTR'
    index_vectr = content.index('VECTR\n')

    # insert the data after VECTR
    for i in range(natom):
        content.insert(index_vectr+1,
                       '{:>5} {:>10.6} {:>10.6} {:>10.6}\n'.format(
                           i+1,
                           scale * np.real(data[interval, point, band, i, 0]),
                           scale * np.real(data[interval, point, band, i, 1]),
                           scale * np.real(data[interval, point, band, i, 2])))
        content.insert(index_vectr+2, '{:>5} 0 0 0 0\n'.format(i+1))
        content.insert(index_vectr+3, '  0 0 0 0 0\n')
        index_vectr += 3

    # find the key word 'VECTT'
    index_vectt = content.index('VECTT\n')

    # insert the data after VECTT
    for i in range(1, natom+1):
        # the number 1 at the end is indicated
        # the middle of arrow place at the atom
        content.insert(index_vectt+1, '{:>5} 0.5 255 0 0 1\n'.format(i))
        index_vectt += 1

    # output the data
    with open('POSCAR_edit.vesta', 'w') as file:
        content = "".join(content)
        file.write(content)


if __name__ == '__main__':
    data = eigenvector_direct_coord('POSCAR', 'band.hdf5')
    replace('POSCAR.vesta', data, 0, 0, 5, 24, scale=30)
