#Damian Rumble, UoE
#20140216

#Run18

#This is my SED plotter - based on 'plotwisesed.py' JH.

#iNPUT CATALOGS - Sorted into various telescopes. X = wavelength (Microns), Y = Flux (mJy)

#####################################################
#import modules

import matplotlib.pyplot as plt
import numpy as np

#####################################################
#constants

mjy = 1e-29
u = 1e-6

########### DATA ###########

#Data - sources from mission specific files which contian flux data for each core on ODD coloumns and error data on EVEN columns.

O = np.loadtxt("MWC297/instruments/optical.txt")
M2 = np.loadtxt("MWC297/instruments/2mass.txt")
S = np.loadtxt("MWC297/instruments/Spitzer.txt")
MSX = np.loadtxt("MWC297/instruments/MSX.txt")
IRAS = np.loadtxt("MWC297/instruments/IRAS.txt")
JCMT = np.loadtxt("MWC297/instruments/JCMT.txt")
SCUBA2 = np.loadtxt("MWC297/instruments/SCUBA2.txt")
IRAM = np.loadtxt("MWC297/instruments/IRAM.txt")
VLA = np.loadtxt("MWC297/instruments/VLA.txt")

########## functions ############

def v(l):
    #simple frequency calculator
    c = 299792458
    return mjy*(c/(l*1E-6))




########## Plotter ############

#set up figure plotter (with subplots)


fig1 = plt.figure(1)
#fig2 = plt.figure(1)
fig1.set_size_inches((7,10))
#fig2.set_size_inches((9,3))
fig1.subplots_adjust(hspace=0.001)
fig1.subplots_adjust(wspace=0.001)
nx=2;ny=5;
fig1.clf()
#fig2.clf()

i = 1
j = 0

for i in np.arange(20):

    if (i % 2): #if index is odd (i.e. flux s selected in the data files as opposed to the flux error)
        plt = fig1.add_subplot(ny,nx,j+1)

        #Plotting SEDS components on to a LogLog Plot complete with Error bars. 

    #Optical
        plt.loglog(O[:,0],O[:,i]*v(O[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(O[:,0],O[:,i]*v(O[:,0]), yerr=O[:,i+1]*v(O[:,0]),fmt='b.', ecolor='blue')

    #2MASS
        plt.loglog(M2[:,0],M2[:,i]*v(M2[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(M2[:,0],M2[:,i]*v(M2[:,0]), yerr=M2[:,i+1]*v(M2[:,0]),fmt='b.', ecolor='blue')

    #Spitzer
        plt.loglog(S[:,0],S[:,i]*v(S[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(S[:,0],S[:,i]*v(S[:,0]),  yerr=S[:,i+1]*v(S[:,0]),fmt='b.', ecolor='blue')

    #MSX
        plt.loglog(MSX[:,0],MSX[:,i]*v(MSX[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(MSX[:,0],MSX[:,i]*v(MSX[:,0]), yerr=MSX[:,i+1]*v(MSX[:,0]),fmt='b.', ecolor='blue')

    #IRAS
        plt.loglog(IRAS[:,0],IRAS[:,i]*v(IRAS[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(IRAS[:,0],IRAS[:,i]*v(IRAS[:,0]),  yerr=IRAS[:,i+1]*v(IRAS[:,0]),fmt='b.', ecolor='blue')

    #JCMT
        plt.loglog(JCMT[:,0],JCMT[:,i]*v(JCMT[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(JCMT[:,0],JCMT[:,i]*v(JCMT[:,0]), yerr=JCMT[:,i+1]*v(JCMT[:,0]),fmt='b.', ecolor='blue')

    #SCUBA2
        plt.loglog(SCUBA2[:,0],SCUBA2[:,i]*v(SCUBA2[:,0]),color='red',marker='x',linestyle='None',markersize=10)
        plt.errorbar(SCUBA2[:,0],SCUBA2[:,i]*v(SCUBA2[:,0]), yerr=SCUBA2[:,i+1]*v(SCUBA2[:,0]),fmt='r.', ecolor='red')
        print SCUBA2[:,i]*v(SCUBA2[:,0])

    #IRAM
        plt.loglog(IRAM[:,0],IRAM[:,i]*v(IRAM[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(IRAM[:,0],IRAM[:,i]*v(IRAM[:,0]),  yerr=IRAM[:,i+1]*v(IRAM[:,0]),fmt='b.', ecolor='blue')

    #VLA
        plt.loglog(VLA[:,0],VLA[:,i]*v(VLA[:,0]),color='black',marker='o',linestyle='None',markersize=2)
        plt.errorbar(VLA[:,0],VLA[:,i]*v(VLA[:,0]),  yerr=VLA[:,i]*v(VLA[:,0]),fmt='b.', ecolor='blue')

        axes=fig1.gca()
        #axes=fig2.gca()

        cores = [1,2,3,4,7,11,12,17,19,30]
        
        for label in axes.get_xticklabels() + axes.get_yticklabels():
            label.set_fontsize('x-small')

        plt.set_xticklabels('',visible='False',size='x-small')

        if (j == 9):
            plt.set_xlabel(r'$\lambda$ (${\mathrm{\mu m}}$)',size='x-small')
            plt.set_xticklabels('',visible='True',size='x-small')
        if (j == 8):
            plt.set_xlabel(r'$\lambda$ (${\mathrm{\mu m}}$)',size='x-small')
            plt.set_xticklabels('',visible='True',size='x-small')

        if (j % 2):
            plt.set_yticklabels('',visible='True',size='x-small')
        else:
            plt.set_ylabel(r'$\nu$ $S_\nu$ (${\mathrm{Wm^{-1}}}$)',size='x-small')

        plt.axis([1e0,2e3,5e-18,5e-10]) #Large Axis
        plt.text(5e2,5e-11,'FW%s'%(cores[j]),size='x-small')

        print 'SED %d plotted'%(cores[j])

        j = j + 1

fig1.savefig('plots/20140327_MWC297_SEDs.eps', format='eps')
