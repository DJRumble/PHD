#Damian Rumble, UoE
#20/04/2015
#alpha.py

#SURFACE plot alpha given beta and temp.

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import cm

def alpha(beta,T):
    ratio = (17**(3+beta)) / (9**(3+beta)) * (np.exp(16.93/T)-1) / (np.exp(31.97/T)-1)
    S850 = 1
    S450 = S850*ratio
    alpha = (np.log10(S450)-np.log10(S850))/0.2769
    return alpha
#axes
T = np.arange(5,45,0.1)
beta = np.arange(0,2.0,0.1)
#grid
beta,T= np.meshgrid(beta, T)
Z = alpha(beta,T)
scale=(beta.min(),beta.max(),beta.min(),beta.max())
#plotting
plt.figure()
im = plt.imshow(Z, interpolation='bilinear', origin='lower',cmap=cm.PRGn, extent=scale)
#colorbar
plt.colorbar(label='alpha')
#lines
plt.axvline(x=1.8,linestyle='dashed',color='red')
plt.axhline(y=0.5,linestyle='dashed',color='red')
#labels
plt.ylabel('Temp. (5K - 45K)')
plt.xlabel('beta')
#contours
levels = (0.0,0.8,1.6)
CS = plt.contour(Z, levels,linewidths=1,colors=('black'),extent=scale)
plt.clabel(CS, levels, inline=1, fmt='%1.1f',fontsize=14)

plt.show()
