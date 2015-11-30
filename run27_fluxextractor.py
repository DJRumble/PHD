#Damian Rumble, UoE
#22/10/2015
#fluxextractor.py

#This code extracts total flux and temp. of a clump, as identified by FellWalker, from a map
#It then calculates:
#- Mass
#- decovolved Radius
#- Column density
#- clump stability


# - distance to the nearest OB star

#All data is written to a file in a suitable format for Latex Tables

#####################################################
#import maths, ploting and astrophysical packages

import numpy as np
import os
import subprocess
import commands
import astropy.io.fits as pyfits
from astropy.io import fits
from astropy import units as u #not used
from astropy.coordinates import SkyCoord  #used
import string

#my modules
import mass_maps2
import jeans
import RaDecDegs
import mass
import ndf2fits
import noise

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

def massF(S,T,kappa,d,mT):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""
    #Equation of Mass in solar masses, per pixel
    if T == 0:
        T = mT
    else:
        M = 1.55*S*(np.exp(17.0/T)-1.0)*((kappa/0.012)**(-1.0))*((d/500.0)**2.0)
    return M

def ColumnDF(M,CD,a,d):
#Convert map of mass (in solar masses) into column density (in H2 cm-2) and Extinction (mag)
    M_x = 1.989E33 #g
    N = 1
    au = 14959787100000 #cm
    m_h = 1.67262178E-24 #g
    mu = 2.8#333  #ratio of H2 to He 5:1

#pixel area
    A = ((((a*d)*au)**2.0)*N) #in cm^2 
    f = M_x/(mu*m_h*A) #column density per cm^2 

    cmd = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,M,f,CD)
    os.system(cmd)
    return

def mass_errF(S,T,dS,dT,mT,mdT): #Mass error
    if T == 0:
        T = mT
    if dT == 0:
        dT = mdT
    S = float(S)
    #v2
    C = (1.55*((kappa/0.012)**(-1.0))*((d/500.0)**2.0))**2.0
    A = ((np.exp(17.0/T)-1.0)*dS)**2.0
    B = (((-17*S)/(T**2.0))*np.exp(17.0/T)*dT)**2.0
    M = C*(A+B)
    dM = np.sqrt(M)
    return dM

########### Bulk Code #############
#deffine constants
kappa = 0.012 #cm^2 g^-1
d = 500 #pc

#Load maps from master file, order: s850,temp,tempr_err,alpha,clumps,clumpcat
maps = np.loadtxt('mapsextractor.txt',dtype='string')

#Open master clump file for writing
SMM = open("SMM/SMM_master.tab","w")
SMM = open("SMM/SMM_master.tab","a")

#File FMT: i,Ra,Dec,s850,M,dM,T,dT,CD,dCD,D,Mj,dMj,M/Mj,dM/Mj,d,N
SMM.write('#i\tRa\tDec\tflux850(Jy)\tMass(Mo)\terrMass\tT(K)\terrT(K)\tCD(1E21_H2cm-2)\tCDerr\tsize(pc)\tdensity(cm-3)\tMj(Mo)\terMj\tM/Mj\terM/Mj\tdistance(pc)\tN\n')


#loop through maps
for i in maps:
    #split INPUT file
    s850,temp,dtemp,alpha,clumps,clumpcat = string.split(i,',')
    prefix = string.split(s850,'.sdf')[0]
    file = string.split(prefix,'/')[1]
    #deffine region name
    name = string.split(clumps,'/')[1]
    region = string.split(name,'_')[0]

    print region#,':',s850,temp,clumps,clumpcat
    
    ### SETUP temp. file directories
    if not os.path.exists('SMM/'+region):
        os.makedirs('SMM/'+region)
    if not os.path.exists('SMM/'+region+'/temp'):
        os.makedirs('SMM/'+region+'/temp')

    ###open maps from FITS files to deffine DIM ###
    dim = PARGET(clumps,'dims','ndftrace')
    a = PARGET(s850,'fpixscale','ndftrace')
    print a[0]

    row = int(dim[1])    #rows of the table from NDFTRACE
    column = int(dim[0]) #columns of the table from NDFTRACE  

    print 'Row: ',row
    print 'Column: ',column

    #CREATE noise map of constant level
    sigma = noise.noise_by_data(s850,'FALSE')
    sigmap = 'SMM/'+region+'/temp/sigmap.sdf'
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%f newhi=%f %s'%(kapdir,s850,sigmap,0,0,sigma,sigma,'> /dev/null')
    os.system(cmd)
    
    #OPEN up the region specific clump list
    data = fits.open(clumpcat)
    DATA = data[1].data

    #LOOP through the clumps in a region
    for j in DATA:
        I = int(j[0])
        print '==========================='
        print 'Clump #'+str(I)+' of '+region
        print '==========================='
        name = region+'-SMM'+str(I)

        ### SETUP temp. file names
        thresh = 'SMM/'+region+'/temp/thresh'+str(I)+'.sdf'
        mask = 'SMM/'+region+'/temp/mask'+str(I)+'.sdf'
        maskbad = 'SMM/'+region+'/temp/maskbad'+str(I)+'.sdf'
        magic = 'SMM/'+region+'/temp/magic'+str(I)+'.sdf'

        ###process 1 - isolate clump i ###
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,clumps,thresh,I,I,0,0,'> /dev/null')
        os.system(cmd)
        ###process 2 - remove 0value data
        cmd = '%s/setmagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
        os.system(cmd)
        ###recreate mask of clump i - mask0
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
        os.system(cmd)
        ###recreate mask of clump i - mask 'bad'
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%d %s'%(kapdir,magic,maskbad,0,0,'bad',1,'> /dev/null')
        os.system(cmd)

        #############
        ####FLUXES###
        #############
        #set up clump file names
        c850 = 'SMM/'+region+'/s850clump'+str(I)+'.sdf'
        cnoise = 'SMM/'+region+'/temp/s850noise'+str(I)+'.sdf'
        #MASK maps
        fmask(mask,s850,c850) #850 flux
        fmask(mask,sigmap,cnoise) #noise mask

        #850#
        cmd = '%s/stats ndf=%s %s'%(kapdir,c850,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        c850i = output
        c850i = round(float(c850i),3)

        print '850um flux = ', c850i,' Jy'

        #################
        ###TEMPERATURE###
        #################
        #masking orginal maps, new saved to temp files
        ctemp = 'SMM/'+region+'/temp/tempr_clump'+str(I)+'.sdf'
        ctemp_err = 'SMM/'+region+'/temp/tempr_clump_err'+str(I)+'.sdf'
        ctempNM = 'SMM/'+region+'/temp/temprNM_clump'+str(I)+'.sdf'
        ctempNM_err = 'SMM/'+region+'/temp/temprNMerr_clump'+str(I)+'.sdf'
        ctempbad = 'SMM/'+region+'/temp/tempr_clump_bad'+str(I)+'.sdf'
        mtemp = 'SMM/'+region+'/tempr_clump'+str(I)+'.sdf'
        mtemp_err = 'SMM/'+region+'/aquila_tempr_clump_err'+str(I)+'.sdf'
        tempNM = 'SMM/'+region+'/temp/temprNM.sdf'
        tempNM_err = 'SMM/'+region+'/temp/temprNM_err.sdf'

        nomagic(temp,tempNM)
        nomagic(dtemp,tempNM_err)

        fmask(maskbad,temp,ctemp) #temp 'bad'
        fmask(maskbad,tempNM,ctempNM) #temp
        fmask(maskbad,dtemp,ctemp_err) #temp_error
        fmask(maskbad,tempNM_err,ctempNM_err) #temp 'bad'

        #replace 'bad' values in temp map with the MEAN temperature of pixels in the clump
        #DATA#
        cmd = '%s/stats ndf=%s %s'%(kapdir,ctemp,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
        status, output = commands.getstatusoutput(cmd)
        temp_mi = float(output)
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
        THtemp = 'SMM/'+region+'/tempr_THclump'+str(I)+'.sdf'
        THtemp_err = 'SMM/'+region+'/tempr_THclump_err'+str(I)+'.sdf'
        cmd = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,ctempNM,THtemp,0.1,1000,temp_mi,0.)
        os.system(cmd)
        cmd = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,ctempNM_err,THtemp_err,0.00001,1000,temperr_mi,0.)
        os.system(cmd)
        #frac error in mean temp
        dTi = temperr_mi/temp_mi

        print 'mean temperature is ',round(temp_mi,1),' pm ',round(temperr_mi,1)#,' pm ',round(temp_sigi,2),' K'

        ############################
        ###Total number of pixels###
        ############################
        cmd = '%s/stats ndf=%s %s'%(kapdir,mask,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        Ni = float(output)

        print 'number of pixels = ',int(Ni)

        ##############################
        #### CALCULATE PROPERTIES ####
        ##############################
        #use NOMAGIC to remove the blanks
        nm850 = 'SMM/'+region+'/temp/s850_NMclump'+str(I)+'.sdf'
        nmtemp = 'SMM/'+region+'/temp/temp_NMclump'+str(I)+'.sdf'
        nmtemp_err = 'SMM/'+region+'/temp/temp_NMclump_err'+str(I)+'.sdf'
        nnoise = 'SMM/'+region+'/temp/s850noiseNM'+str(I)+'.sdf'
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
        nm850FITS = 'SMM/'+region+'/temp/s850_NMclump'+str(I)+'.fits'
        cnoise = 'SMM/'+region+'/temp/s850noiseNM'+str(I)+'.fits'
        nmtempFITS = 'SMM/'+region+'/temp/temp_NMclump'+str(I)+'.fits'
        nmtempFITS_err = 'SMM/'+region+'/temp/temp_NMclump_err'+str(I)+'.fits'
        #new files#
        massFITS = 'SMM/'+region+'/Mass'+str(I)+'.fits'
        masserrFITS = 'SMM/'+region+'/masserr'+str(I)+'.fits'

        #copy existing clump file to act as a base on to which calculations will be writen.
        copy(cnoise,massFITS)
        copy(cnoise,masserrFITS)
 
        #open map data
        I850,D850 = openimage(nm850FITS) #open flux data
        Itemp,Dtemp = openimage(nmtempFITS) #open temp data
        Itemperr,Dtemperr = openimage(nmtempFITS_err) #open temp error
        Inoise,Dnoise = openimage(cnoise) #open noise
        #ICD,CD = openimage(CDFITS)  #CD data in a table
        Imass,mass = openimage(massFITS)  #mass data in a table
        Imasserr,masserr = openimage(masserrFITS)  #masserr data in a table

        #run property calculations per pixel via various functions. 
        #Where there is no temp data set values to 0.
        k = 0
        l = 0
        #print D850[0][k][l]
        #print Dtemp[k][l]
        for k in range(0, row):
            for l in range(0, column):
                S = float(D850[k][l])
                T = float(Dtemp[k][l])
                dS = float(Dnoise[k][l])
                dT = float(Dtemperr[k][l])
                if S == 0:
                    mass[k][l] = 0
                    masserr[k][l] = 0
                if S != 0:
                    mass[k][l] = massF(S,T,kappa,d,temp_mi)
                    masserr[k][l] = mass_errF(S,T,dS,dT,temp_mi,temperr_mi)
        Imasserr.flush()
        Imass.flush()

        ######
        #MASS#
        ######
        #convert back to sdf
        ndf2fits.fits2ndf(massFITS)
        ndf2fits.fits2ndf(masserrFITS)
        mass = 'SMM/'+region+'/Mass'+str(I)+'.sdf'
        masserr ='SMM/'+region+'/masserr'+str(I)+'.sdf'

        #Extract TOTAL clump mass as number
        cmd = '%s/stats ndf=%s %s'%(kapdir,mass,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        mi = float(output)
        del output
       
        #calculate frac error on the mass and find its mean value, on which the mean mass error is based. 
        fracmass = 'SMM/'+region+'/temp/fracmass_clump'+str(I)+'.sdf'
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
        CD ='SMM/'+region+'/ColumnD'+str(I)+'.sdf'

        #calculate CD map from individual clumps
        ColumnDF(mass,CD,a[0],d)

        #Run Stats and extract TOTAL value of CD
        cmd = '%s/stats ndf=%s %s'%(kapdir,CD,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'MAXIMUM','stats') 
        status, output = commands.getstatusoutput(cmd)
        Avi = float(output)/(0.9E21) #H2 per cm2
        CDi = float(output)/1E21
        #calculate CD error from mass fractional error
        dCDi = fMi*CDi
        print 'Max CD = ',round(CDi,2),'pm',round(dCDi,2),' x1E21 H2 cm-2'
        print 'max extinction = ',round(Avi,2)

        size = (j[5]+j[6])/2.
        size_cm = (size*3.)*500*1.49597877E13 #convert units to cm
        size_pc = size_cm/3.08567758E18
        print 'size = ',round(size_pc,2),' psc'
        density = (CDi*1E21)/size_cm
        print 'Density = ',round(density,-4), ' cm-3'

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





        ############
        #WRITE FILE#
        ############
        SMM.write(str(name)+'\t'+str(j[1])+'\t'+str(j[2])+'\t'+str(c850i))
        SMM.write(str(mi)+'\t'+str(dMi)+'\t'+str(temp_mi)+'\t'+str(temperr_mi))
        SMM.write('\t'+str(CDi)+'\t'+str(dCDi)+'\t'+str(size_pc)+'\t'+str(density))
        SMM.write('\t'+str(Ni)+'\n')

    break
