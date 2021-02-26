import numpy as np

data = np.fromfile("wtf.raw", dtype=np.int16)
print(data.size)
import matplotlib.pyplot as plt

plt.plot(data)
plt.show()