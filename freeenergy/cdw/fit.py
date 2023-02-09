import numpy as np
import matplotlib.pyplot as plt

data_x = np.random.normal(0, 1, 10)
data_y = np.random.normal(0, 1, 10)

plt.scatter(data_x, data_y)
fit = np.polyfit(data_x, data_y, 3)
t = np.poly1d(fit)
x = np.linspace(0, 0.5, 100)
plt.plot(t(x))

plt.show()
