# Damian Rumble, UoE
# 20131211
# Convert decimal degrees to starlink grid coords

import os
import re
import astropy.io.fits as pyfits
import numpy as np
import RaDecDegs

stardir = '/star/bin/'

#params = 'outparams.fit'
#refmap = 'SerpensMWC297_20130918_232648_s850_DJR_cal.sdf'



def coord(refmap,params):
    #takes input reference map (.sdf) and parameter file as produced by FellWalker) (.fit)
    data = pyfits.getdata(params)

    array = []

    i = 0
    for i in range(len(data)):

        ra = RaDecDegs.timeRa(data[i][1],'TRUE')
        dec = RaDecDegs.timeDec(data[i][2],'TRUE')
    
        coords = "%10.7f %s"%(ra, dec)
        cmd = stardir+"kappa/wcstran %s posin=\"\'%s\'\" framein=SKY frameout=PIXEL > /dev/null"%(refmap,coords)
        os.system(cmd)
        posout = (os.popen(stardir+"kappa/parget posout wcstran").readlines())
    #print posout[0]
        (pixx,pixy,junk) = re.split("\s+",posout[0],2)
    #print pixx, pixy

        #print 'Ra '+str(ra)+' and Dec '+ str(dec) +' are coordinates '+ str(pixx)+' and ' + str(pixy)
        
        array.append([pixx,pixy])
        i = i + 1
    return array

if __name__ == "__main__":
    coord(refmap,params)


