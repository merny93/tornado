import matplotlib.pyplot as plt
import numpy as np
import spify 

a1 = [[4.9500,4.9454,4.9455,4.941],[0.0004,0.0004,0.0004,0.003]]
a2 = [[5.8,5.76,5.77,5.827],[0.2,0.02,0.03,0.003]]
c2 = [[3.0,3.19,3.17,3.180],[0.2,0.02,0.02,0.002]]

x = [0,25,50,75,100]

figs, axs = plt.subplots(3)
axs[0].scatter(x[:4], a1[0], marker = '.')
axs[0].errorbar(x[:4], a1[0], yerr = a1[1], linestyle = 'None')

axs[1].scatter(x[1:], a2[0], marker = '.')
axs[1].errorbar(x[1:], a2[0], yerr = a2[1], linestyle = 'None')
axs[2].scatter(x[1:], c2[0], marker = '.')
axs[2].errorbar(x[1:], c2[0], yerr = c2[1], linestyle = 'None')
axs[0].set_ylabel(r'$a_{fcc}$')
axs[1].set_ylabel(r'$a_{tetra}$')
axs[2].set_ylabel(r'$c_{tetra}$')
plt.xlabel('Tin %')
for ax in axs:
    ax.set_xlim(-1,101)
plt.savefig('tinleadalloy.png')