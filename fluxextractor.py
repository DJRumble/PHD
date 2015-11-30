#Damian Rumble, UoE
#22/11/2013
#run14_fluxextractor.py

#This code extracts total flux of a clump, as identified by FellWalker, and the temperature for the same region and stores it in a file for calculation as a mass.

#Flux typically has more pixels than temp. Where there is overlap, bad temp pixels are labelled with the mean temp. of the remaining group. Where no temp. data exists, pixels are labelled with a constant temperature of 15K 

#Jeans Mass is calculated in parrallel 

#All data is written to a file in a suitable format for Latex Tables

#####################################################
#import maths, ploting and astrophysical packages

import numpy as np
import os
import subprocess
import commands
import astropy.io.fits as pyfits
from astropy import units as u #not used
from astropy.coordinates import SkyCoord  #used
import string

#my modules
import mass_maps2
import jeans
import RaDecDegs
import mass
import ndf2fits

#####################################################
#constants

kappa_8 = 0.012 #opacity
kappa_4 = 0.048 #opacity
d = 500 #distance (pcs)
errd = 50 #distance error
#each pixel has an area in SI units and each apature contains N pixels
a = 2.999997
#au = 149597871000 #m
#pi = 3.14159265359
#psc = 3.08567758E+16 #m
#m_h2 = 3.34745E-24 #g
#mu = 2.333333  #ratio of H2 to He (5:1)

########### set CSH Commands directories #############

kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

########### Functions #############

def nomagic(infile,outfile):
    cmd = '%s/nomagic in=%s out=%s repval=0 QUIET'%(kapdir,infile,outfile)
    os.system(cmd)

def copy(infile,outfile):
    cmd = 'cp %s %s'%(infile,outfile)
    os.system(cmd)

def openimage(file):
    image = pyfits.open(str(file), mode = 'update')
    data = image[0].data
    return image,data

def fmask(mask,inmap,outmap):
    mul1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,mask,inmap,outmap)
    os.system(mul1)
    #print 'masked'
    return outmap

def PARGET(file,parameter,kappa):
    #Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file.
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)
    a = np.loadtxt('parameter.txt')
    os.remove('parameter.txt')
    return a

def massF(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""
    #Equation of Mass in solar masses, per pixel
    if T == 0:
        T = 15
    M = 1.55*S*(np.exp(17.0/T)-1.0)*((kappa/0.012)**(-1.0))*((d/500.0)**2.0)
    return M

def ColumnDF(M,a,d):
    #constants
    au = 1.49597871E13 #cm
    #au = 149597871000 #m
    M_x = 1.989E33 #g
    m_h = 1.67262178E-24 #g
    mu = 2.333333  #ratio of H to He (5:1)
    A = ((((a*d)*au)**2.0)) #in cm
    m = M*M_x #A = 1 density in g per pixel
    n = m/(mu*m_h*A) #column density per cm^2 
    #N = m/A #density in g per cm^2
    return n

def mass_errF(S,T,dS,dT): #Mass error
    if T == 0:
        T = 15.
    if dT == 0:
        dT = 2.
    S = float(S)
    #v2
    C = (1.55*((kappa/0.012)**(-1.0))*((d/500.0)**2.0))**2.0
    A = ((np.exp(17.0/T)-1.0)*dS)**2.0
    B = (((-17*S)/(T**2.0))*np.exp(17.0/T)*dT)**2.0
    M = C*(A+B)
    dM = np.sqrt(M)
    return dM

def CD_errF(CD,mass,mass_err): #CD error
    if mass == 0:
        errCD = 0
    else:
        frac = mass_err / mass
        errCD = frac * CD
    return errCD**2.




########### Bulk Code #############
#deffine constants
kappa = 0.012 #cm^2 g^-1
d = 500 #pc

#clumps = 'output/Aquila_extS2_850-4-30am_clumps5.sdf'
clumps = 'output/S2-FF/Aquila-FF_extS2_850-4-30am_clumps5.sdf'
clumpsH70 = 'output/Aquila_850_clumps3_H70.sdf'

#deffine input maps
#SDF
s450 = 'input/S2-ff/Aquila_extS2_450-4-30am_fullfreefree.sdf'
s850 = 'input/S2-ff/Aquila_extS2_850-4-30am_fullfreefree.sdf'
temp = 'input/S2-ff/Aquila_extS2-4-30am_fullfreefreetemperatureNM.sdf'
temp_err = 'input/S2-ff/Aquila_extS2-4-30am_fullfreefreetemp_errorNMth.sdf'
noise  = 'input/Aquila_850noise.sdf'
H70 = 'input/aquilaM2-070.sdf'
YSO = 'output/surface/Aquila850_YSOsurface_r300_units.sdf'
YSOpoint = 'output/surface/YSOmap850.sdf'

###open maps from FITS files to deffine DIM ###
#deffine input maps - use clump map as it is a different size to the main maps (resulting from ndfsections) - all masks are based in this map.
dim = PARGET(clumps,'dims','ndftrace')

row = int(dim[1])    #rows of the table from NDFTRACE
column = int(dim[0]) #columns of the table from NDFTRACE  
print 'row = ' + str(row) 
print 'column = ' + str(column)  

a = 3 #pixel size

### Deffine Catalog files and headers ###
#open files for writing
SMM = open("SMM/SMM-Aq-ff.tab","w")
SMM = open("SMM/SMM-Aq-ff.tab","a")
SMM.write('#i\tRa\tDec\tflux450(Jy)\tflux850(Jy)\tfluxH70(MJy/Sr)\tMass(Mo)\terrMass\tT(K)\terrT(K)\tCD(1E22_H2cm-2)\tCDerr\tYSO\terrYSO\tsumYSO\tMj(Mo)\terMj\tM/Mj\terM/Mj\tdistance(pc)\tN\n')

#Open latex readable file 1
SMM_LTX = open("SMM/SMM-Aq-ff_LTX.txt","w")
SMM_LTX = open("SMM/SMM-Aq-ff_LTX.txt","a")
SMM_LTX.write('index\t&\tCoords\t&\t\tflux450\t&\tflux850\t&\tfluxH70\t&\tdistance\n')
SMM_LTX.write('\t&\t(J2000)\t&\t\t(Jy)\t&\t(Jy)\t&\t(MJy/Sr)\t&\t(pc)\n')

#Open Latex readable file 2
SMM_LTX2 = open("SMM/SMMprop-Aq-ff_LTX.txt","w")
SMM_LTX2 = open("SMM/SMMprop-Aq-ff_LTX.txt","a")
SMM_LTX2.write('index\t&\tflux850\t&\tMass\t&\tTemp.\t&\tColumn Density\t&\tYSOdensity\t&\tYSOs\t&\tMj\t&\tM/Mj\n')
SMM_LTX2.write('index\t&\t(Jy)\t&\t(Mo)\t&\t(K)\t&\t(1E22 H2cm-2)\t&\t(YSO per pc-2)\t&\t(per clump)\t&\t(Mo)\t&\t\n') 

### REGRIDDING ###
#For 450 clumps, the 450 map needs to be regridded onto 850 scale to match the clumps (also 850 scale). Then the total flux devided by 2.25 to regain 450 levels. This requires collapsing first.
d2_450 = 'temp/s450_2d.sdf'
d2_850 = 'temp/s850_2d.sdf'
align = 'temp/Aquila_extS2_450-4-30am_align_850.sdf'

cmd = '%s/collapse in=%s out=%s axis=%s %s'%(kapdir,s450,d2_450,'3','> /dev/null')
os.system(cmd)
cmd = '%s/collapse in=%s out=%s axis=%s %s'%(kapdir,s850,d2_850,'3','> /dev/null')
os.system(cmd)
cmd = '%s/wcsalign in=%s out=%s ref=%s method=%s conserve=%s %s'%(kapdir,d2_450,align,d2_850,'nearest','TRUE','accept')
os.system(cmd)

### Extract FW clumps
#params = 'output/Aquila_s850_SerpSouth_clump_cat.tab'
#params = 'output/Aquila_s850_W40_clump_cat.tab'
#params = 'output/AquilaFF_s850_SerpSouth_clump_cat.tab'
#params = 'output/AquilaFF_s850_W40_clump_cat.tab'	
params = 'output/S2-FF/AquilaFF_clump_cat.tab'

coord = np.loadtxt(params)

flux450 = []
flux850 = []

i = 0
l = len(coord[:,0])

#Start clump loop - cylce each clump, mask S2 data, compute properties and write to file. 
for i in range(l):
    I = int(coord[i][0])
    print '##############################'
    print 'Starting processing of clump#'+str(I)
    print '##############################'
    if i == 200:
        print 'Not a real clump'
    else:

    #Set up temp. file names
        thresh = 'temp/thresh'+str(I)+'.sdf'
        threshH70 = 'temp/threshH70'+str(I)+'.sdf'
        magic = 'temp/magic'+str(I)+'.sdf'
        mask = 'temp/mask'+str(I)+'.sdf'
        magicH70 = 'temp/magicH70'+str(I)+'.sdf'
        maskH70 = 'temp/maskH70'+str(I)+'.sdf'
        maskbad = 'temp/maskbad'+str(I)+'.sdf'
        mtemp = 'temp/meantemp'+str(I)+'.sdf'
        ctempbad = 'temp/temp_clump_bad'+str(I)+'.sdf'

        ###process 1 - isolate clump i ###
        #note that i on the H70 clump map is the integral of i*ratio of pixel areas
        factor = 1.13569 #(3.19706^2)/(3^2)
        IH70 = int(I*factor)
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,clumps,thresh,I,I,0,0,'> /dev/null')
        os.system(cmd)
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,clumpsH70,threshH70,IH70,IH70,0,0,'> /dev/null')
        os.system(cmd)
    
        ###process 2 - remove 0value data
        cmd = '%s/setmagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
        os.system(cmd)
        cmd = '%s/setmagic in=%s out=%s repval=%d %s'%(kapdir,threshH70,magicH70,0,'> /dev/null')
        os.system(cmd)

        ###recreate mask of clump i - mask0
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
        os.system(cmd)
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magicH70,maskH70,0,0,0,1,'> /dev/null')
        os.system(cmd)
 
        ###recreate mask of clump i - mask 'bad'
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%d %s'%(kapdir,magic,maskbad,0,0,'bad',1,'> /dev/null')
        os.system(cmd)

        #############
        ####FLUXES###
        #############
        #set up clump file names
        c450 = 'clumps/450/aquila_s450_clump'+str(I)+'.sdf'
        c850 = 'clumps/850/aquila_s850_clump'+str(I)+'.sdf'
        cH70 =  'clumps/H70/aquila_H70_clump'+str(I)+'.sdf'
        cnoise = 'temp/Aquila_850noise'+str(I)+'.sdf'
        #masking orginal maps, new saved to temp files
        fmask(mask,align,c450) #450 flux
        fmask(mask,s850,c850) #850 flux
        fmask(mask,noise,cnoise) #noise mask

        #Extract total fluxes
        #450#
        cmd = '%s/stats ndf=%s %s'%(kapdir,c450,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        c450i = output 
        c450i = round(float(c450i),2)
        #850#
        cmd = '%s/stats ndf=%s %s'%(kapdir,c850,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        c850i = output
        c850i = round(float(c850i),3)
        #H70#
        fmask(maskH70,H70,cH70) #70um flux
        #Extract total fluxes
        cmd = '%s/stats ndf=%s %s'%(kapdir,cH70,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        cH70i = float(output) 
        
        #print results#
        print '850um flux = ', c850i,' Jy'
        print '70um flux = ', cH70i,' MJy/Sr' 

        ############################
        ###Total number of pixels###
        ############################
        cmd = '%s/stats ndf=%s %s'%(kapdir,mask,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        Ni = float(output)

        print 'number of pixels = ',int(Ni)

        #################
        ###TEMPERATURE###
        #################
        #masking orginal maps, new saved to temp files
        ctemp = 'temp/aquila_temp_clump'+str(I)+'.sdf'
        ctemp_err = 'temp/aquila_temp_clump_err'+str(I)+'.sdf'
        ctemp_var = 'temp/aquila_temp_clump_var'+str(I)+'.sdf'
        mtemp = 'clumps/temperature/aquila_temp_clump'+str(I)+'.sdf'
        mtemp_err = 'clumps/temperature/aquila_temp_clump_err'+str(I)+'.sdf'
        #fmask(maskbad,temp,ctempbad) #temp 'bad'
        fmask(maskbad,temp,ctemp) #temp
        fmask(maskbad,temp_err,ctemp_err) #temp_error
        
        #replace 'bad' values in temp map with the MEAN temperature of pixels in the clump
        #DATA#
        cmd = '%s/stats ndf=%s %s'%(kapdir,ctempbad,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
        status, output = commands.getstatusoutput(cmd)
        temp_mi = float(output)
        #Run Stats and extract STDV value of the temp
        #cmd = '%s/parget parname=%s applic=%s'%(kapdir,'sigma','stats') 
        #status, sigma = commands.getstatusoutput(cmd)
        #temp_sigi = float(sigma)/(np.sqrt(Ni))
        #ERROR#
        #cmd = '%s/pow in=%s power=%s out=%s'%(kapdir,ctemp_err,2.,ctemp_var) 
        #os.system(cmd)
        cmd = '%s/stats ndf=%s %s'%(kapdir,ctemp_err,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
        status, output = commands.getstatusoutput(cmd)
        temperr_mi = (float(output))

        if temp_mi == 0:
            temp_mi = 15.0
        if temperr_mi == 0:
            temperr_mi = 2.0

        #replace '0' values with 'mean' values
        THtemp = 'temp/aquila_temp_THclump'+str(I)+'.sdf'
        THtemp_err = 'temp/aquila_temp_THclump_err'+str(I)+'.sdf'
        cmd = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,ctemp,THtemp,0.1,1000,temp_mi,0.)
        os.system(cmd)
        cmd = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,ctemp_err,THtemp_err,0.00001,1000,temperr_mi,0.)
        os.system(cmd)

        #frac error in mean temp
        dTi = temperr_mi/temp_mi

        print 'mean temperature is ',round(temp_mi,1),' pm ',round(temperr_mi,1)#,' pm ',round(temp_sigi,2),' K'

        ##############################
        #### CALCULATE PROPERTIES ####
        ##############################
        #use NOMAGIC to remove the blanks
        c850 = 'temp/aquila_s850_NMclump'+str(I)+'.fits'
        nm850 = 'temp/aquila_s850_NMclump'+str(I)+'.sdf'
        nmtemp = 'temp/aquila_temp_NMclump'+str(I)+'.sdf'
        nmtemp_err = 'temp/aquila_temp_NMclump_err'+str(I)+'.sdf'
        nnoise = 'temp/Aquila_850noiseNM'+str(I)+'.sdf'
        nomagic(c850,nm850)
        nomagic(cnoise,nnoise)
        nomagic(THtemp,nmtemp)
        nomagic(THtemp_err,nmtemp_err)

        #from individual clump maps
        ndf2fits.ndf2fits(nm850)
        ndf2fits.ndf2fits(nnoise)
        ndf2fits.ndf2fits(nmtemp)
        ndf2fits.ndf2fits(nmtemp_err)      
        #deffine FITS files
        #existing files#
        cnoise = 'temp/Aquila_850noiseNM'+str(I)+'.fits'
        nmtempFITS = 'temp/aquila_temp_NMclump'+str(I)+'.fits'
        nmtempFITS_err = 'temp/aquila_temp_NMclump_err'+str(I)+'.fits'
        #new files#
        CDFITS = 'clumps/CD/Aquila-ColumnD'+str(I)+'.fits'
        massFITS = 'clumps/mass/Aquila-Mass'+str(I)+'.fits'
        masserrFITS = 'clumps/mass/Aquila-masserr'+str(I)+'.fits'
        #copy existing clump file to act as a base on to which calculations will be writen.
        copy(c850,CDFITS)
        copy(c850,massFITS)
        copy(c850,masserrFITS)
 
        #open map data
        I850,D850 = openimage(c850) #open flux data
        Itemp,Dtemp = openimage(nmtempFITS) #open temp data
        Itemperr,Dtemperr = openimage(nmtempFITS_err) #open temp error
        Inoise,Dnoise = openimage(cnoise) #open noise
        ICD,CD = openimage(CDFITS)  #CD data in a table
        Imass,mass = openimage(massFITS)  #mass data in a table
        Imasserr,masserr = openimage(masserrFITS)  #masserr data in a table

        #run property calculations per pixel via various functions. 
        #Where there is no temp data set values to 0.
        j = 0
        k = 0
        for j in range(0, row):
            for k in range(0, column):
                S = float(D850[j][k])
                T = float(Dtemp[j][k]) #temp_mi #
                dS = float(Dnoise[j][k])
                dT = float(Dtemperr[j][k])
                if S == 0:
                    mass[j][k] = 0
                    masserr[j][k] = 0
                    CD[j][k] = 0
                    #errCD[j][k] = 0
                if S != 0:
                    mass[j][k] = massF(S,T,kappa,d)
                    masserr[j][k] = mass_errF(S,T,dS,dT)
                    CD[j][k] = ColumnDF(float(mass[j][k]),a,d)

        ICD.flush()
        Imasserr.flush()
        Imass.flush()

        #convert back to sdf
        ndf2fits.fits2ndf(CDFITS)
        ndf2fits.fits2ndf(massFITS)
        ndf2fits.fits2ndf(masserrFITS)
        CD = 'clumps/CD/Aquila-ColumnD'+str(I)+'.sdf'
        mass = 'clumps/mass/Aquila-Mass'+str(I)+'.sdf'
        masserr = 'clumps/mass/Aquila-masserr'+str(I)+'.sdf'
        #mCD = 'clumps/CD/aquila_CD+m_clump'+str(I)+'.sdf'

        ######
        #MASS#
        ######
        #Extract TOTAL clump mass as number
        cmd = '%s/stats ndf=%s %s'%(kapdir,mass,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        mi = float(output)
        del output
        #calculate frac error on the mass and find its mean value, on which the mean mass error is based. 
        fracmass = 'temp/fracmass_clump'+str(I)+'.sdf'
        cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,masserr,mass,fracmass)
        os.system(cmd)
        cmd = '%s/stats ndf=%s %s'%(kapdir,fracmass,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
        status, output = commands.getstatusoutput(cmd)
        dMi = np.sqrt(float(output))
        fMi = (float(output))
        dMi = fMi*mi
        print 'mass = ',round(mi,2),'pm',round(dMi,2),' Mdot: ', round(fMi*100.,2),'%'

        ####
        #CD#
        ####
        #Run Stats and extract TOTAL value of CD
        cmd = '%s/stats ndf=%s %s'%(kapdir,CD,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'MAXIMUM','stats') 
        status, output = commands.getstatusoutput(cmd)
        CDi = float(output)/1E22
        #calculate CD error from mass fractional error
        dCDi = fMi*CDi
        print 'Max CD = ',round(CDi,2),'pm',round(dCDi,2),' 1E22 H2 cm-2'

        #####################
        #YSO surface density#
        #####################
        #files
        cYSO = 'clumps/YSOsurface/aquila-YSOdensity'+str(I)+'.sdf'
        pYSOi = 'clumps/YSOsurface/aquila-YSO'+str(I)+'.sdf'
        #masking
        fmask(mask,YSO,cYSO)
        fmask(mask,YSOpoint,pYSOi)
        #Run Stats and extract MEAN value of the surface density
        cmd = '%s/stats ndf=%s %s'%(kapdir,cYSO,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
        status, output = commands.getstatusoutput(cmd)
        YSO_mi = float(output)
        #Run Stats and extract STDV value of the surface density
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'sigma','stats') 
        status, sigma = commands.getstatusoutput(cmd)
        YSO_sigi = float(sigma)/(np.sqrt(Ni))
        #extract the total numbner of YSO from the clump
        cmd = '%s/stats ndf=%s %s'%(kapdir,pYSOi,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        nYSOi = float(output)
        print 'Mean YSO density is ',round(YSO_mi,1),'pm',round(YSO_sigi,1),' in YSOs per pc^2'
        print 'Number of YSOs is ',int(nYSOi)

        ############
        #jeans mass#
        ############
        #Run Stats and extract NUMBER of pixels
        cmd = '%s/stats ndf=%s %s'%(kapdir,magic,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'numgood','stats') 
        status, output = commands.getstatusoutput(cmd)
        Ni = float(output)

        #run module for jeans mass and error
        Ji = jeans.sadavoy(temp_mi,d,Ni)
        err_Ji = Ji * dTi
        
        #run calculation of jeans stability and error
        if float(Ji) == 0:
            equ = 0
        else:
            equ = float(mi) / float(Ji)

        #Error in mass ratio:
        frac = ((float(dMi)/float(mi))**2.)+((err_Ji/float(Ji))**2.)
        err_equ = equ * ((frac)**0.5)

        print 'Jeans mass = ',round(Ji,2),'pm',round(err_Ji,2),' Mdot'
        print 'Jeans mass ratio = ',round(equ,2),'pm',round(err_equ,2)

        ############################
        #Distance from OS1a - in pc#
        ############################
        OS1ra = 277.865945833
        OS1dec = -2.08968333333
        rai = float(coord[i][1])
        deci = float(coord[i][2])
      
        c = SkyCoord(ra=rai*u.degree, dec=deci*u.degree, distance=d*u.pc, frame='icrs')
        o = SkyCoord(ra=OS1ra*u.degree, dec=OS1dec*u.degree, distance=d*u.pc, frame='icrs')
        r = o.separation_3d(c)
        R = string.split(str(r),' pc')[0]
        R = round(float(R),2)
        print 'distance to OS1a = ',R

        ###########################
        #Round data to sig figures#
        ###########################
        c850i = round(c850i,2)
        cH70i = round(cH70i,-6)
        mi = round(mi,2)
        dMi = round(dMi,2)
        temp_mi = round(temp_mi,1)
        temperr_mi = round(temperr_mi,1)
        CDi = round(CDi,3)
        dCDi = round(dCDi,3)
        YSO_mi = round(YSO_mi,0)
        YSO_sigi = round(YSO_sigi,0)
        Ji = round(Ji,1)
        err_Ji = round(err_Ji,1)
        equ = round(equ,1)
        err_equ = round(err_equ,1)

    #RAW data file
        SMM.write(str(I)+'\t')
        SMM.write(str(round(rai,4))+'\t')
        SMM.write(str(round(deci,4))+'\t')
        SMM.write(str(c450i)+'\t')
        SMM.write(str(c850i)+'\t')
        SMM.write(str(cH70i)+'\t')
        SMM.write(str(mi)+'\t')
        SMM.write(str(dMi)+'\t')
        SMM.write(str(temp_mi)+'\t')
        SMM.write(str(temperr_mi)+'\t')
        SMM.write(str(CDi)+'\t')
        SMM.write(str(dCDi)+'\t')
        SMM.write(str(YSO_mi)+'\t')
        SMM.write(str(YSO_sigi)+'\t')
        SMM.write(str(Ji)+'\t')
        SMM.write(str(err_Ji)+'\t')
        SMM.write(str(equ)+'\t')
        SMM.write(str(err_equ)+'\t')
        SMM.write(str(R)+'\n')

    #Convert Ra/Dec from degrees to Time format.
        c = SkyCoord(ra=rai*u.degree, dec=deci*u.degree)
        hms = c.to_string('hmsdms')
        print hms

    #LaTex data file
        SMM_LTX.write('Aq-SMM'+str(I)+'\t&\t')
        SMM_LTX.write('JCMTLSG J'+str(hms)+'\t&\t')
        SMM_LTX.write(str(c450i)+'\t&\t')
        SMM_LTX.write(str(c850i)+'\t&\t')
        SMM_LTX.write(str(cH70i)+'\t&\t')
        SMM_LTX.write(str(R)+'aaa\n')

    #LaTex data file
        SMM_LTX2.write('Aq-SMM'+str(I)+'\t&\t')
        SMM_LTX2.write(str(c850i)+'\t&\t')
        SMM_LTX2.write(str(mi)+'$\pm$'+str(dMi)+'\t&\t')
        SMM_LTX2.write(str(temp_mi)+'$\pm$'+str(temperr_mi)+'\t&\t')
        SMM_LTX2.write(str(CDi)+'$\pm$'+str(dCDi)+'\t&\t')
        SMM_LTX2.write(str(YSO_mi)+'$\pm$'+str(YSO_sigi)+'\t&\t')
        SMM_LTX2.write(str(Ji)+'$\pm$'+str(err_Ji)+'\t&\t')
        SMM_LTX2.write(str(equ)+'$\pm$'+str(err_equ)+'aaa\n')
    #increment
    i = i + 1
    
print 'complete'
