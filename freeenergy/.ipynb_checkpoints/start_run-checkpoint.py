from free import get_dos
from free import do_integral
from free import fermienergy
from free import plot_dos
from free import get_mu
from free import thermointegral
from free import E, dos, T_num
import numpy as np
import matplotlib.pyplot as plt


def plot():
    phi_data = []
    T_data = []
    phi_output = ""

    [dos_data, E, dos] = get_dos('dos_input')
    K0 = do_integral(E, dos)
    minus_K0 = -1*K0
    K2 = do_integral(E, minus_K0)
    [Ef, Ef_index, n0] = fermienergy(1, 3E-29, K0)

    plot_dos(dos_data)

    Tmin = 0
    Tmax = 0
    [mu_data, mu_index_data] = get_mu(n0, K0, Tmin, Tmax, Ef, Ef_index)

    for i in range(T_num+1):
        T = Tmin + (Tmax - Tmin)*i/T_num
        phi = thermointegral(K2, T, mu_data[i], mu_index_data[i], Ef)
        phi_data[i] = [T, phi + n0*mu_data[i]]
        phi_output = phi_output+str(T)+"\t"+str(phi)+"\n"
        T_data[i] = T
