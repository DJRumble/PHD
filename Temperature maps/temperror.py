#Damian Rumble, UoE
#28/05/2013
#run8.py

#This is a simple code to numerically calculate the error in temperature and create .txt file to be used in TOPCAT

import numpy as np
import pyfits
import sys

##### DEFFINE Parrameters ##########

beta = 1.8     #beta value held constant for calculations
lookup = np.zeros((99501, 20), float)
T = 5.0

#### Sandbox for analytical error solution #######
'''
beta = 1.8 
T = 100
dT = 15
S_r = 10.3685195532

X = dT/T
Y = np.exp(16.93/T)
Z = np.exp(31.97/T)
V = np.exp(-15.04/T)
U = np.exp(-31.97/T)

A = (17**(3+beta)) / (9**(3+beta))

#solution by derrivaitve of Expotnets and division
dS_r1 = S_r*(((((-16.93*(X**2)*Y)/(Y-1))**2)+(((-31.97*(X**2)*Z)/(Z-1))**2))**0.5)

#solution by quotient rule
dS_r2 = A*dT*(((15.04*Y)+(16.93*V)-31.97)/((T**2)*(Z+U-2)))
'''

##### Main Program #################

dT = 0.25 

file = open('temp25.txt', 'w')
#write header
file.write('#T(K) dT(K) S_r dS_r % +dT -dT +dS_r -dS_r \n')
#For range of temp. of  from 5.0 going up in increments of 0.01, specific ratios are calculated.
for i in range(0, 99501):
    lookup[i][0] = float(T)
    #creates temperature-ratio lookup table
    T += 0.01 

    #set temp error (%) 
    lookup[i][1] = dT*lookup[i][0]

    #creates plus and minus error margins
    lookup[i][4] = (1+dT)*lookup[i][0]
    lookup[i][5] = (1-dT)*lookup[i][0]

    #creates ratio from temp
    lookup[i][2] = (17**(3.0+beta)) / (9**(3.0+beta)) * (np.exp(16.93/lookup[i][0])-1) / (np.exp(31.97/lookup[i][0])-1)
    #creates ratio error from temp error
    lookup[i][6] = (17**(3.0+beta)) / (9**(3.0+beta)) * (np.exp(16.93/lookup[i][4])-1) / (np.exp(31.97/lookup[i][4])-1)
    lookup[i][7] = (17**(3.0+beta)) / (9**(3.0+beta)) * (np.exp(16.93/lookup[i][5])-1) / (np.exp(31.97/lookup[i][5])-1)
    #Mean error in ratio
    lookup[i][3] = ((lookup[i][6]-lookup[i][2]) + (lookup[i][2]-lookup[i][7]))/2.0

##### analytical solution  #################

    X = lookup[i][1]/lookup[i][0]
    Y = np.exp(16.93/lookup[i][0])
    Z = np.exp(31.97/lookup[i][0])
    V = np.exp(-15.04/lookup[i][0])
    U = np.exp(-31.97/lookup[i][0])

    A = (17**(3+beta)) / (9**(3+beta))

#solution by derrivaitve of Expotnets and division
    lookup[i][8] = lookup[i][2]*(((((-16.93*(X**2)*Y)/(Y-1))**2)+(((-31.97*(X**2)*Z)/(Z-1))**2))**0.5)

#solution by quotient rule
    lookup[i][9] = A* lookup[i][1]*(((15.04*Y)+(16.93*V)-31.97)/((lookup[i][0]**2)*(Z+U-2)))




##### write to text file #################
    file.write(str(lookup[i][0])) #temperature/K
    file.write(' ')
    file.write(str(lookup[i][1])) #temp error
    file.write(' ')
    file.write(str(lookup[i][2])) #Fluxratio
    file.write(' ')
    #file.write(str(lookup[i][3])) #mean flux error
    #file.write(' ')
    #file.write(str(100*(lookup[i][3]/lookup[i][2]))) #mean flux error
    #file.write(' ')
    #file.write(str(lookup[i][8])) #Analytical solution (1)
    #file.write(' ')
    #file.write(str(100*(lookup[i][8]/lookup[i][2]))) #mean flux error
    #file.write(' ')
    file.write(str(lookup[i][9])) #Analytical solution (2-quotient)
    file.write(' ')
    file.write(str(100*(lookup[i][9]/lookup[i][2]))) #mean flux error
    file.write(' ')
    file.write(str(lookup[i][4])) #+sigT
    file.write(' ')
    file.write(str(lookup[i][5])) #-sigT
    file.write(' ')
    file.write(str(lookup[i][6])) #+sigR 
    file.write(' ')
    file.write(str(lookup[i][7])) #-sigR
 

    file.write('\n')


