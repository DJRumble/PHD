import numpy as np

h = 6.62607004E-34
c = 2.998E8
v = 3.52697E11
k = 1.3806488E-23
T = 20
pc = 2.0856E16
kappa = 0.0035 
Jy = 1E-26
S = 0.2612 #half times two   #0.0171 #small  #0.367 #Large
d = 250

MJ = 1.898E27

B = (2.*h*(v**3.)/(c**2.))*(1/(np.exp((h*v)/(k*T))-1))

M = (S*Jy*(d*pc)**2.)/(kappa*B)

m = M/MJ

print 'Mass of disc = ',m,' in Mj' 
