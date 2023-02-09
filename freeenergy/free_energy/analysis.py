import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from free import create_plot
from free import load
import ipdb


if __name__ == "__main__":
    Ef = 4.039503297421883

    data = load("data_dos.csv")
    energy = data['energy']
    dos = data['dos']

    fig, ax = create_plot('density of state')
    ax.plot(energy-Ef, dos)
    ax.fill_between(energy-Ef, dos, where=energy < Ef, color='gray')
    plt.show()
