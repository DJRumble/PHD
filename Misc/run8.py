#Damian Rumble, UoE
#24/06/2013
#run8.py

#A sand box for plotting results for errors in temp.

import numpy as np
from pylab import *


##### DEFFINE Constants ##########

#dx is % error in temp. a constant, derived at T=13.5K, likewise beta

fdT1 = 0.03 #3%
fdT2 = 0.10 #10%
fdT3 = 0.25 #25%
fdT4 = 0.33 #33%

beta = 1.8

beta = 0.0

for i in range(200):
    beta = beta + 0.01
    A = (17.0**(3.0+beta)) / (9.0**(3.0+beta))




A = (17.0**(3.0+beta)) / (9.0**(3.0+beta))

##### DEFFINE Variables ##########

#x axis = temperature
#y axis = ratio error

T = arange(5.0, 80.0, 0.01)
dT = fdT1*T

##### DEFFINE component functions ##########

def Y(T):
    Y = np.exp(16.93/T) 
    #print Y
    return Y

def Z(T):
    Z = np.exp(31.97/T)
    #print Z
    return Z

def V(T):
    V = np.exp(-15.04/T)
    return V

def U(T):
    U = np.exp(-31.97/T)
    return U

##### DEFFINE mainfunctions ##########

#z = Z(T)
#y = Y(T)
#print Z(T)
#print Y(T)


def S_R(T,dT):
    # flux ratio (T)
    S_R = A * ((Z(T)-1)/(Y(T)-1))
    plot(dT,S_R)
    return

def del_S_R(T,dT):
    # error in flux ratio (T,dT)
    del_S_R  = A * dT * (((15.04 * Y(T))+(16.93 * V(T))-31.97)/((T*T) * (Z(T)+U(T)-2.0)))
    S_R = A * ((Z(T)-1)/(Y(T)-1))
    plot(dT,del_S_R)
    return




##### Plotting ##########


#S_R(T,dT)
del_S_R(T,dT) 

##### Axis labels etc ##########

xlabel('Temp. error')
ylabel('Flux ratio error')
title('dT_dS_03')
grid(True)
savefig("temp.pdf")
show()
