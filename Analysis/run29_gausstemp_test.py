

#DJR UoE
#30/10/2015 - run29.py
#
#code is designed to calculate the gassian wieghted temp of a core. 

###########  Import modules ######################################

import os
import numpy as np
import string
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
import copy
import commands
from astropy import units as u
from astropy.coordinates import SkyCoord
import parget
import re

########### set CSH Commands directories #############

kapdir = '/stardev/bin/kappa'
convdir = '/stardev/bin/convert'

###################################################################
#Start of bulk script
###################################################################

def exp(x,y,FWHM,out,null):
    ln = 2.354820045 #2SQRT(2ln2)
    p = 3
    sig = 2*(((FWHM/p)/ln)**2.0)
    exp = "'exp(-(((xa+0.5"+x+")**2)+((xb+0.5"+y+")**2))/"+str(sig)+")'"

    print exp

    #norm = 1/(sig*np.sqrt(2*np.pi))
    #print norm
    
    cmd = '%s/maths exp=%s out=%s like=%s'%(kapdir,exp,out,null)
    os.system(cmd)
    #cmd = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,out,norm,Nout)
    #os.system(cmd) 
    return 

def guasmooth(x,y,fwhm,pix):
    pixel = "'Gausstests/pixel.sdf'"

    #Implant point sources into null maps
    section = "\"'%i,%i'\""%(int(x),int(y))
    chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,null,pixel,section,n)
    os.system(chpix1)  
    sig = fwhm/(2.*(2.*np.log(2))**0.5)
    smooth = sig/3
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,pixel,smooth,out)
    os.system(cmd)
    return 

def fmask(maskmap,inmap,outmap):
    #create mask and mask input map
    maskbad = 'MSK.sdf'  
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%d %s'%(kapdir,maskmap,maskbad,0,0,'bad',1,'> /dev/null')
    os.system(cmd)
    mul1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,maskbad,inmap,outmap)
    os.system(mul1)
    #print 'masked'
    return outmap

#Produce YSO map
maps = np.loadtxt('maps_gausstemp.txt',dtype='string')

#Open master results file for writing
table = open("GBSprotostar_temperatures_test.tab","w")
table = open("GBSprotostar_temperatures_test.tab","a")
table.write('#region\tra(deg)\tdec(deg)\tGra(deg)\tGdec(deg)\ttemp(K)\ttemperr(K)\talphaE\tTbolE\n')

#LOAD protostars
pos = np.loadtxt('GBS_YSO_test.tab',dtype='string')

#loop through maps
for i in maps:
    #split INPUT file - #declare map origin
    s850,temp,temperr,raO,decO,d = string.split(i,',')
    prefix = string.split(s850,'.sdf')[0]
    file = string.split(prefix,'/')[1]
    #deffine region name
    name = string.split(s850,'/')[1]
    region = string.split(name,'_')[0]
    print region
    
    #Produce null map
    null = 'null.sdf'
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,null,0,0,0,0)
    os.system(thresh1)

    #deffine FWHM
    core = 0.05 #pc
    FWHM = (((core*3.08567758E16)/1.49597871E11)/float(d)) #arcsecs
    print 'FWHM is '+str(round(FWHM,1))

    #declare degrees per pixel per wavelength.
    pix850 = [8.3333333E-4,-8.33333333E-4]
    #pix850 = [9E-3,-8.3333333333E-4]

    n = 1
    k = 0
    l = 0

    #deffine limits on region search area
    limitra,limitdec = 1,1

    ralo = float(raO)-float(limitra)
    rahi = float(raO)+float(limitra)
    declo = float(decO)-float(limitdec)
    dechi = float(decO)+float(limitdec)
    print ralo, rahi, declo, dechi
    print raO,decO

    for j in pos:
        #limit catalogue to the specific region
        if (float(j[2]) > ralo) & (float(j[2]) < rahi) & (float(j[3]) > declo) & (float(j[3]) < dechi):
            if float(j[4]) > -0.3:
                RaH = float(j[2])/15
                DecD = float(j[3])
                coords = "%10.7f %s"%(RaH, DecD)

                cmd = "%s/wcstran ndf=%s posin=\"\'%s\'\" framein=SKY frameout=PIXEL > /dev/null"%(kapdir,s850,coords)
                print cmd
                os.system(cmd)
                posout = (os.popen(kapdir+"/parget posout wcstran").readlines())
                #posout = parget.PARGET(s850,'posout','wcstran')
                #print posout
                (x,y,junk) = re.split("\s+",posout[0],2)
                #break

                #calculate displacement of source from the origin
                delx,dely = (float(raO)-float(j[2])),(float(decO)-float(j[3]))
                #list pixel displacements
                xold,yold = (delx/pix850[0]),(dely/pix850[1])

                print 'OLD',round(xold,1),round(yold,0)

                if x < 0:
                    X = '+'+str(abs(int(round(float(x),1))))
                if y < 0:
                    Y = '+'+str(abs(int(round(float(y),1))))
                if x > 0:
                    X = '-'+str(int(round(float(x),1)))
                if y > 0:
                    Y = '-'+str(int(round(float(y),1)))

                X = str(int(round(float(x),1)))
                Y = str(int(round(float(y),1)))

                print 'NEW',X,Y

                #Produce gaussian 
                out = 'Gausstests/out.sdf'#%(j[1])
                Gmask = 'Gmask.sdf'
                #exp(X,Y,FWHM,out,null)
                guasmooth(X,Y,FWHM,pix850[0])
                #mask gaussian
                fmask(temp,out,Gmask)

                #SUM gauss component
                cmd = '%s/stats ndf=%s %s'%(kapdir,Gmask,'> /dev/null')
                os.system(cmd)
                cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
                status, output = commands.getstatusoutput(cmd)
                gw = float(output)

                #Record gaussian position 
                cmd = '%s/stats ndf=%s %s'%(kapdir,out,'> /dev/null')
                os.system(cmd)
                
                status, output = commands.getstatusoutput(cmd)
                ysoloc = str(output)
                ysoloc = string.split(ysoloc,',')
                YSOloc = ysoloc[0]+ysoloc[1]
                c = SkyCoord('%s'%(YSOloc), unit=(u.hourangle, u.deg))
                
                difX = (float(c.ra.degree)-float(raO))/float(pix850[0])
                difY = (float(c.dec.degree)-float(decO))/float(pix850[1])

                table.write(str(out)+'\t'+str(j[2])+'\t'+str(j[3])+'\t'+str(c.ra.degree)+'\t'+str(c.dec.degree)+'\t'+str(round(difX,1))+'\t'+str(round(difY,1))+'\n')
                break

    break


