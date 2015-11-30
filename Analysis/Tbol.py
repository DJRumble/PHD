#Damian Rumble, UoE
#20140401

#Tbol

#This code/module calculates Tbol across an SEDS - with errors via numerics

#iNPUT CATALOGS - Sorted into various telescopes. X = wavelength (Microns), Y = Flux (mJy)

#####################################################
#import modules

import matplotlib.pyplot as plt
import numpy as np
import random as RD

import trapizum

####################################################
#constants

mjy = 1e-29
bol = 1.25E-11
d = 250*(3.085678E16)
L = 3.839E26
pi = 3.141592

########### DATA ###########

#Data - data for specific YSO (bolometric and submm). l(0),S_v(1),dS_v(2)


#bolometric = np.loadtxt("bolometric.txt")
#SubMM = np.loadtxt("SubMM.txt")
#enoch = np.loadtxt("enoch_test.txt")

########## functions ############

def v(l):
    #simple frequency calculator
    c = 299792458
    return (c/l)

T = open("bolometric_temperatures.txt","w")
T = open("bolometric_temperatures.txt","a")

T.write('#YSO\t&T_bol(K)\t&L_bol(Lo)\t&\tL_r(%)\n')

def TBOL(nu,S,Snu):
    #calculate integrals for I0 and I1
    SUM_Snu = -trapizum.trap(nu,Snu)
    SUM_S = -trapizum.trap(nu,S)
    V_m = SUM_Snu/SUM_S
    TBOL = V_m*bol
    return TBOL

def LBOL(d,nu,nu_smm,S,S_smm):
    L = 3.839E26

    #calculate integrals for I0 and I1
    SUM_S = -trapizum.trap(nu,S)

    #SubMM integral
    SMM = -trapizum.trap(nu_smm,S_smm)

    #luminosity calc.
    L_bol = (4*np.pi*(d**(2.))*SUM_S)/L
    L_SMM = (4*np.pi*(d**(2.))*SMM)/L

    L_R = (L_SMM/L_bol)*100

    Lumin = [L_bol,L_SMM,L_R]

    return Lumin

j = 0
k = 0

F = 0
f = 0
Fnu = 0 

YSO = [1,2,3,4,7,11,12,17,19,30]

i = 0

for i in range(len(YSO)):

    ####Bulk calculations####

    print '=============='
    print 'YSO #', YSO[i]
    print '=============='
    print 'Loading files:'

    file_enoch = 'MWC297/instruments/enoch_test.txt'
    file_dunham = 'MWC297/instruments/dunham_test.txt'
    file_bol =  'MWC297/YSO/bol/FW'+"%i"%(YSO[i])+'_bol.txt'
    file_smm =  'MWC297/YSO/smm/FW'+"%i"%(YSO[i])+'_smm.txt'
    file_smm_dunham =  'MWC297/instruments/dunham_smm_test.txt'

    file_bol = np.loadtxt(str(file_bol))
    file_smm = np.loadtxt(str(file_smm))
    file_enoch = np.loadtxt(str(file_enoch))
    file_dunham = np.loadtxt(str(file_dunham))
    file_smm_dunham = np.loadtxt(str(file_smm_dunham))

    print 'Processing catalog:'

    nu = v(file_bol[:,0]*1E-6)  
    S = file_bol[:,1]*mjy 
    Snu = S*nu

    nu_smm = v(file_smm[:,0]*1E-6)  
    S_smm = file_smm[:,1]*mjy 

    #Calculate TBol
    Tbol = TBOL(nu,S,Snu)

    #Calculate LSMM/LBol             
    Lbol = LBOL(d,nu,nu_smm,S,S_smm)

    ####Error calculation####
    
    Tbol_dist = []
    Lbol_dist = []
    Lsmm_dist = []
    Lr_dist = []

    for j in range(100000):
        S = (RD.normalvariate(file_bol[:,1],file_bol[:,2]))*mjy
    
        Tbol_RNG = TBOL(nu,S,Snu)
        Lbol_RNG = LBOL(d,nu,nu_smm,S,S_smm)
     
        Tbol_dist.append(Tbol_RNG)
        Lbol_dist.append(Lbol_RNG[0])
        Lsmm_dist.append(Lbol_RNG[1])
        Lr_dist.append(Lbol_RNG[2])

    errorTbol = np.std(Tbol_dist)
    errorLbol = np.std(Lbol_dist)
    errorLr = np.std(Lr_dist)


    ###Printing results###

    print 'Calculating TBOL:'
    print 'Tbol = ',"%1.f"%Tbol,' (',"%1.f"%(errorTbol), ') K'

    print 'Calculating LBOL:'  
    print 'Lbol = ',"%.2f"%(Lbol[0]),' (',"%.2f"%(errorLbol), ') L_dot'
    #print 'LSMM = ',Lbol[1],,'pm',"%.4f"%(errorLbol), ' L_dot'
    #print 'ratio (bol/smm) = ',"%.2f"%(Lbol[0]/Lbol[1])
    print 'ratio (%) = ',"%.2f"%(Lbol[2]),' (',"%.2f"%(errorLr),')'

    T.write("%i"%(YSO[i])+'\t&\t')
    T.write("%1.f"%Tbol+' '+"%1.f"%(errorTbol)+'\t&\t')
    T.write("%.4f"%(Lbol[0])+' '+"%.2f"%(errorLbol)+'\t&\t')
    T.write("%.2f"%(Lbol[2])+' '+"%.2f"%(errorLr)+'\n')
