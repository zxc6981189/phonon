import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ipdb

T_num = 300
kb = 1.380658E-23  # Bolzmann in Joule/Kelvin
ec = 1.60217733E-19  # 1 eV in Joul
kB = kb/ec

volume = 557.39E-30
ne = 192


def __init__():
    """
    initial setting for matplotlib
    """
    # global font style setting
    plt.rcParams.update({'font.family': 'sans-serif'})
    plt.rcParams.update({'font.serif': 'Helvetica'})
    plt.rc('font', weight='bold')
    plt.rc('font', size=18)


def get_dos(filename='dos_input'):
    """read dos from file"""
    # read without header and use whitespace
    dos = pd.read_csv(filename, header=None, delim_whitespace=True)
    dos = dos.astype(float)
    dos = dos.to_numpy()

    # vasp use state/eV
    dos_data = [dos[:, 0], dos[:, 1]/volume]
    E = dos[:, 0]
    dos = dos[:, 1]/volume
    return [dos_data, E, dos]


def do_integral(x, y):
    """
    do the integration in the interval of x to x+1
    args:
        x(np.array): list of data along x axis
        y(np.array): list of data of f(x)
    """
    result = np.zeros(x.shape[0])
    length_x = x.shape[0]
    length_y = y.shape[0]

    # the length of  data x and y should be the same
    if length_x == length_y:
        for i in range(1, length_x):
            result[i] = result[i-1] + 0.5*(y[i]+y[i-1])*(x[i]-x[i-1])

    return result


def find_index(value, index_guess, array):
    """
    look for the nearest index in the array
    args:
        value(number): desire value
        index_guess(int): guess desire value
        array(np.array): data array

    return:
        int: the desire index in the array
    """

    if array.shape[0] <= index_guess:
        index_guess = array.shape[0]-1

    if array[index_guess] < value:
        while array[index_guess] < value:
            index_guess += 1
            if index_guess >= array.shape[0]:
                break
        index = index_guess-1

    else:
        while array[index_guess] > value:
            index_guess -= 1
            if index_guess < 0:
                break
        index = index_guess

    return index


def fermienergy(E, dos, ne, volume, K0):
    """
    find out the fermienergy in the system
    args:
        E(np.array): energy list
        dos(np.array): dos list
        ne(int): number of electron in the primitive cell
        volume(float): primitive cell volume
        K0(np.array): accumulate dos list

    return:
        fermienergy, index of fermienergy in the nearest energy data list
    """
    n0 = ne/volume

    Ef_index = find_index(n0, K0.shape[0], K0)

    Ef = E[Ef_index] + 2*(n0-K0[Ef_index])/(dos[Ef_index+1]+dos[Ef_index])

    return [Ef, Ef_index, n0]


def interval(E, mu, T, i, K):
    """
    args:
        E(np.array): energy data list
        mu(number): guess mu
        T(number): temperature
        i(int): interval of temperature index
        K(np.array): accumulate dos

    return:
        electron density in guess mu
    """
    x1 = (E[i]-mu)/(kB*T)
    x2 = (E[i+1]-mu)/(kB*T)
    b = (K[i+1]-K[i])/(E[i+1]-E[i])
    It = (b*(E[i+1]-mu)*np.exp(x2)+b*(E[i]-mu)-K[i]) / \
        (np.exp(x2)+1)-kB*T*b*np.log(np.exp(x2)+1)
    It = It - (b*(E[i]-mu)*np.exp(x1)+b*(E[i]-mu)-K[i]) / \
        (np.exp(x1)+1)+kB*T*b*np.log(np.exp(x1)+1)
    return It


def thermointegral(K, T, mu, mu_index, E, Ef):
    """
    args:
        K(np.array): accumulate dos
        T(number): temperature
        mu(number): chemical potential
        mu_index(int): chemical potential index in the energy list
        E(np.array): energy list
        Ef(number): fermienergy

    return:
        grand potential phi
    """
    if T == 0:
        quantity = K[mu_index]+(K[mu_index+1]-K[mu_index]) * \
            (Ef-E[mu_index])/(E[mu_index+1]-E[mu_index])
    else:
        quantity = 0
        qmax = np.round(64*kB*T/(E[mu_index+2]-E[mu_index]))+1
        for q in range(int(-qmax), int(qmax)+1):
            quantity = quantity + interval(E, mu, T, mu_index+q, K)
    return quantity


def get_mu(n0, K0, Tmin, Tmax, E, Ef, Ef_index):
    mu_data = np.zeros(T_num+1, dtype=np.float64)
    mu_index_data = np.zeros(T_num+1, dtype=np.int64)
    mu_start = Ef*0.9
    mu_old = Ef
    mu_new = Ef
    mu_start_index = find_index(mu_start, Ef_index, E)
    mu_old_index = find_index(mu_old, Ef_index, E)
    mu_new_index = find_index(mu_new, Ef_index, E)
    error_start = n0 - \
        thermointegral(K0, Tmin, mu_start, mu_start_index, E, Ef)
    for i in range(T_num+1):
        T = Tmin + (Tmax - Tmin)*i/T_num
        error_old = n0 - thermointegral(K0, T, mu_old, mu_old_index, E, Ef)
        while (np.abs(error_old/n0) > 1E-7):
            mu_new = mu_old - error_old * \
                (mu_old-mu_start)/(error_old-error_start)
            mu_new_index = find_index(mu_new, mu_old_index, E)
            error_start = error_old
            error_old = n0 - thermointegral(K0, T, mu_new, mu_new_index, E, Ef)
            mu_start = mu_old
            mu_old = mu_new
            mu_old_index = mu_new_index
        mu_data[i] = mu_new
        mu_index_data[i] = mu_new_index
    return [mu_data, mu_index_data]


def create_plot(title, nrows=1, ncols=1):
    """
    This function make the plot to be consistent.
    """
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 10))
    fig.suptitle(title)
    return fig, ax


def save(data, filename):
    """
    args:
        data(np.array): the data desire to save
        filename(str): file name to save the data
        header(bol): (default: False) save data withou header
        index(bol): (default: False) save data withou row index
    """
    data = pd.DataFrame(data)
    data.to_csv(filename, index=False)


def load(filename):
    """
    load data from file
    args:
        filename(str): name of data file
    return:
        data(dataframe): data read from csv file
    """
    data = pd.read_csv(filename)
    return data


def main():
    free = np.zeros(T_num+1)
    temp = np.zeros(T_num+1)

    [dos_data, E, dos] = get_dos('DOSCAR')
    K0 = do_integral(E, dos)
    minus_K0 = -1*K0
    K2 = do_integral(E, minus_K0)
    [Ef, Ef_index, n0] = fermienergy(E, dos, ne, volume, K0)

    Tmin = 0
    Tmax = 600
    [mu_data, mu_index_data] = get_mu(n0, K0, Tmin, Tmax, E, Ef, Ef_index)

    for i in range(T_num+1):
        T = Tmin + (Tmax - Tmin)*i/T_num
        phi = thermointegral(K2, T, mu_data[i], mu_index_data[i], E,  Ef)

        free[i] = phi + n0*mu_data[i]
        temp[i] = T

    print('Ef: ', Ef)
    data_free = {'temp': temp, 'phi': free}
    data_dos = {'energy': E, 'dos': dos}
    save(data_free, 'data_free.csv')
    save(data_dos, 'data_dos.csv')


if __name__ == '__main__':
    main()
