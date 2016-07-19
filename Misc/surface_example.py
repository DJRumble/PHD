import numpy as np
from pylab import *

# the function that I'm going to plot
def z_func(T,dT):
    #ranges
    #T = arange(0.1, 80.0, 1)
    #dT = fdT1*T
    beta = 1.8
    A = (17.0**(3.0+beta)) / (9.0**(3.0+beta))
    Y = np.exp(16.93/T)
    Z = np.exp(31.97/T)
    V = np.exp(-15.04/T)
    U = np.exp(-31.97/T)
    del_S_R = A * dT * (((15.04 * Y)+(16.93 * V)-31.97)/((T*T) * (Z+U-2.0)))
    return del_S_R
 
def z_mass(s,t):
    #ranges
    #t1 = np.arange(0.1, 50.0, 0.1)
    #s1 = np.arange(0.0005, 0.25, 0.0005)
    A = 5.03517447875E28 #m2
    kappa = 0.00191
    return 0.19*s*A*(np.exp(17/t)-1)*((kappa/0.01)**(-1))


fdT1 = 0.03
t1 = np.arange(0.1, 50.0, 1)
s1 = np.arange(0.0005, 0.25, 0.0005)

X,Y = meshgrid(s1,t1) # grid of point
Z = z_func(X, Y) # evaluation of the function on the grid

im = imshow(Z) # drawing the function
# adding the Contour lines with labels
cset = contour(Z,arange(0,1.3,0.3),linewidths=2,cmap=cm.Set1)
clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
colorbar(im) # adding the colobar on the right
# latex fashion title
grid(True)
#title('Surface plot of error in Flux ratio')
xlabel('Temp. (K)')
ylabel('Temp. error (K)')
show()
