import numpy as np
import matplotlib.pyplot as plt
deg = 1
data = np.loadtxt("calib_data.txt", delimiter=" ")
data = np.roll(data, 1, axis=1)
# print(type(data))
data_lin = data[5:20, :]
plt.plot(data[:,0], data[:,1])
plt.plot(data_lin[:,0], data_lin[:,1])
p = np.polyfit(data_lin[:,0], data_lin[:,1],deg)
x_fine = np.linspace(data_lin[0,0], data_lin[-1,0])
d_pred = np.polyval(p, data_lin[:,0])
plt.plot(x_fine, np.polyval(p, x_fine))
plt.show()

chi = np.sum(((d_pred - data_lin[:,1])/0.05)**2)
print(chi/(d_pred.size - deg - 1))
print("polynomial is", p)
plt.scatter(np.arange(d_pred.size), d_pred - data_lin[:,1] )
plt.show()