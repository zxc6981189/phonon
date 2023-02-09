#! /Users/liangying/Documents/Code/FreeEnergy/.venv/bin/python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import ipdb


def hdf5_parse(filename='band.hdf5'):
    """
    parse the hdf5 file

    return:
        data: hdf5 data
    """
    data = h5py.File(filename)
    return data


def create_plot():
    # global font style setting
    plt.rcParams.update({'font.family': 'serif'})
    plt.rcParams.update({'font.serif': 'Times New Roman'})
    plt.rcParams.update({'font.size': 22})

    # fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 10))
    fig.tight_layout()
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.3,
                        hspace=0.3)
    return fig, axes


def plot(data, color='blue', filename='phonon_dispersion.png'):
    """
    plot phonon dispersion and save as png file

    args:
        data(hdf5_parse): the object return from hdf5_parse function
        color(string): band line color default in blue
        filename(str): output figure filename
    return:
        output the saved phonon dispersion figure
    """

    # create the plot
    fig, axis = create_plot()

    # number of intervals between two q points
    # number of points in the interval
    # total number of band in the system
    intervals, points, nband = data['frequency'].shape

    # eigenvalue aka phonon frequency
    frequency = data['frequency']

    # the actual q points distance from the frist point of the path
    distance = data['distance']

    # the label of each interval
    label_interval = data['label'][:].astype('str').tolist()

    # change the label_interval into a 1d list without the duplicate labels
    labels = [label_interval[0][0]]

    # the x position of the label
    # this is the list for plotting
    label_position = [0]
    for i in range(len(label_interval)):
        labels.append(label_interval[i][-1])
        label_position.append(distance[i][-1])

    # figure q points range
    x_min = np.min(distance[:][4])
    x_max = np.max(distance[:][4])

    # figure frequency range
    # define by hand
    y_min = -4
    y_max = 10

    # zero frequency line
    axis.plot([x_min, x_max], [0, 0], '--', c='red')
    axis.set_xlim(x_min, x_max)

    for index_line in range(1, intervals):
        axis.plot([label_position[index_line], label_position[index_line]],
                  [y_min, y_max], '-', c='grey')

    # plot phonon bands
    for index_intervals in range(intervals):
        for index_nband in range(nband):
            x = distance[index_intervals]
            y = frequency[index_intervals, :, index_nband]

            # color in bule (default)
            axis.plot(x, y, color)

    # set xticts with q vector label
    axis.set_xticks(label_position, labels=labels)

    # the range in y axis is dependent on different case
    axis.set_ylim(y_min, y_max)
    axis.set_xlim(x_min, x_max)

    fig.savefig(filename, transparent=True)


if __name__ == '__main__':
    # outuput the figure
    data = hdf5_parse()
    plot(data, sys.argv[1])
