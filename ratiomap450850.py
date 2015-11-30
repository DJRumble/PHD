#################################################################################
###########################   RATIOMAP450850.py   ###############################
#################################################################################

#Damian Rumble, UoE
#20141130

#This code takes SCUBA-2 data at 450um and 850um and prodcues maps of flux ratio.
#These are then used to create temperature maps at a fixed beta. 

#This code is structured as follows:

#A) Free Parameters
#B) Fixed parameters 
#1) Structuring of directories
#2) 4 component dual beam convolution
#3) Alignment 
#4) Threshing - calculates noise level in script from module 'noise.py'
#5) Making the ratio maps
#6) Making uncertianty maps
#7) Producing temperture maps from look up table
#8) Making uncertainty maps of temperature
#9) Cleaning up

#################################################################################
########### import modules  #############
#Standard modules
import numpy as np
#import astropy.io.fits as pyfits
import pyfits
import sys
import os
import string
import mpfit
import copy

#My modules
import noise #My noise calculating module


#################################################################################
########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'
#################################################################################
########### Functions
def PARGET(file,parameter,kappa):
    #Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file.
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)

    a = np.loadtxt('parameter.txt')
    os.remove('parameter.txt')
    return a

def FWHM(value,LOGIC):
    #LOGIC is 'TRUE' or 'FALSE'.
    #LOGIC == TRUE. Transform FWHM into Sigma
    #LOGIC == FALSE. Transform Sigma into FWHM
    if LOGIC == 'TRUE':
        out = value/(2.*(2.*np.log(2))**0.5)
    elif LOGIC == 'FALSE':
        out = value * (2.*(2.*np.log(2))**0.5)
    else:
        print 'Please enter a LOGIC value as "TRUE" or "False"'
    return out

def zoominratio(i,j,min, max, count,ratio,lookup):
    mid = (min+max)/2  #sets middle
    count += 1         #increments flag
    if(count < 14):
        #if value in top half of search, search again with top half of section of array
        if(ratio[i][j] > lookup[mid][1]):
            zoominratio(i,j,int(mid)-1, max, count,ratio,lookup)
        #if in bottom half, call bottom half of array again
        elif(ratio[i][j] < lookup[mid][1]): 
            zoominratio(i,j,min, int(mid)+1, count,ratio,lookup)
    #once zoomed into small enough section of array, find best fit
    else: 
        result = 1                                                             
        difference = 100
        for k in range(min, max):
            #check if closest match or not
            if (np.abs(lookup[k][1] - ratio[i][j]) < difference):
                result = k #set array index if smallest so far
                difference = np.abs(lookup[k][1] - ratio[i][j]) #update smallest difference 
        ratio[i][j] = lookup[result][0] 
    #output correct temperature once best fit found
    return

def temperature(sdf,fits,var,beta,LOGIC):
    #Calculate dim numbers from ndftrace 
    dim = PARGET(sdf,'dims','ndftrace')
    #number of dimensions, this item contains to numbers, typically 571 columns 536 row
##### DEFFINE Parrameters ##########
    image = pyfits.open(str(fits), mode = 'update') #opens FITS file created from NDF
    ratio = image[0].data  #Ratio flux data in a table
    if (LOGIC == 'TRUE'):
        variance =  image[1].data  #opens variance array directly from the ndf array
    elif (LOGIC == 'FALSE'):
        imvariance = pyfits.open(str(var), mode = 'update') #opens Variance FITS file 
        variance =  imvariance[0].data #variance data in a table
    row = int(dim[1])    #rows of the table from NDFTRACE
    column = int(dim[0]) #columns of the table from NDFTRACE   
    lookup = np.zeros((99501, 3), float)         
    upperratio = np.zeros((row, column), float)
    lowerratio = np.zeros((row, column), float)
    T = 5.0 #minimum temperature 
##### Main Program ################
    print 'Beta = ' + str(beta)
    print 'row = ' + str(row)
    print 'column =' + str(column)
    print 'image =' + str(image)
#Look up table is created for range of temp. from 5.0K going up in increments of 0.01K.
    for i in range(0, 99501):
        lookup[i][0] = float(T)
        T += 0.01
        lookup[i][1] = (17**(3+beta)) / (9**(3+beta)) * (np.exp(16.93/lookup[i][0])-1) / (np.exp(31.97/lookup[i][0])-1) #Temperature Equation#
#compare table of ratios with real data and use this to build new temp maps - based on varience
    for i in range(0, row):
        for j in range(0, column):
        #if variance flag set, find upper and lower limits based on ratio variance
            if (ratio[i][j] > 0 and LOGIC == "TRUE"):
                upperratio[i][j] = ratio[i][j] + variance[i][j]
                lowerratio[i][j] = ratio[i][j] - variance[i][j]
    for i in range(0, row):
        for j in range(0, column):
            if ratio[i][j] > 0:
                zoominratio(i,j,0, 99500, 0,ratio,lookup) #call functions that find best temperature fit for ratios

#output new values (temperature from ratio) back to FITS file and save
    image.flush() 
    image.close()           

def temp_error(Sratio,tempCUT,ferror,Temp_error,Temp_percent):
    #Naming some files. All in MATHs directory which is subsequently deleted to save on space
    th17 = output_dir+'/math/17.sdf'
    th32 = output_dir+'/math/32.sdf'
    th_15 = output_dir+'/math/neg15.sdf'
    th_32 = output_dir+'/math/neg32.sdf'
    y = output_dir+'/math/y.sdf'
    Y = output_dir+'/math/Y.sdf'
    z = output_dir+'/math/z.sdf'
    Z = output_dir+'/math/Z.sdf'
    v = output_dir+'/math/v.sdf'
    V = output_dir+'/math/V.sdf'
    u = output_dir+'/math/u.sdf'
    U = output_dir+'/math/U.sdf'
    sqT = output_dir+'/math/sqT.sdf'
    Z_1 = output_dir+'/math/z_1.sdf'
    Y_1 = output_dir+'/math/Y_1.sdf'
    Z_Y = output_dir+'/math/Z_Y.sdf'
    ZU = output_dir+'/math/ZU.sdf'
    ZU_2 = output_dir+'/math/ZU_2.sdf'
    sq_T_ZU2 = output_dir+'/math/sqrt2.sdf'
    numerator = output_dir+'/math/numerator.sdf'
    Y15 = output_dir+'/math/15Y.sdf'
    V17 = output_dir+'/math/17V.sdf'
    YV = output_dir+'/math/YV.sdf'
    dinominator = output_dir+'/math/dinomintor.sdf'
    fraction = output_dir+'/math/fraction.sdf'
    Frac_err = output_dir+'/math/frac_err.sdf'
    
    # I need to create number maps for the division for each constant involved here.
    print '-Threshing'
    cmd1 = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,Sratio,th17,0,0,16.93,16.93)
    cmd2 = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,Sratio,th32,0,0,31.97,31.97)
    cmd3 = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,Sratio,th_15,0,0,-15.04,-15.04)
    cmd4 = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,Sratio,th_32,0,0,-31.97,-31.07)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    print '-Functions'
    cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,th17,tempCUT,y)
    os.system(cmd)
    cmd = '%s/expon in=%s out=%s base=Natural'%(kapdir,y,Y)
    os.system(cmd)
    cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,th32,tempCUT,z)
    os.system(cmd)
    cmd = '%s/expon in=%s out=%s base=Natural'%(kapdir,z,Z)
    os.system(cmd)
    cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,th_15,tempCUT,v)
    os.system(cmd)
    cmd = '%s/expon in=%s out=%s base=Natural'%(kapdir,v,V)
    os.system(cmd)
    cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,th_32,tempCUT,u)
    os.system(cmd)
    cmd = '%s/expon in=%s out=%s base=Natural'%(kapdir,u,U)
    os.system(cmd)
    #numerator (top) - Arithmatic
    print '-Numerator'
    cmd1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,tempCUT,tempCUT,sqT)
    cmd2 = '%s/csub in=%s scalar=1 out=%s'%(kapdir,Z,Z_1)
    cmd3 = '%s/csub in=%s scalar=1 out=%s'%(kapdir,Y,Y_1)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,Z_1,Y_1,Z_Y)
    os.system(cmd)
    cmd = '%s/add in1=%s in2=%s out=%s'%(kapdir,Z,U,ZU)
    os.system(cmd)
    cmd = '%s/csub in=%s scalar=2 out=%s'%(kapdir,ZU,ZU_2)
    os.system(cmd)
    cmd = '%s/mult in1=%s in2=%s out=%s'%(kapdir,ZU_2,sqT,sq_T_ZU2)
    os.system(cmd)
    cmd = '%s/mult in1=%s in2=%s out=%s'%(kapdir,sq_T_ZU2,Z_Y,numerator)
    os.system(cmd)
    print '-Dinominator'
    cmd1 = '%s/cmult in=%s scalar=15.04 out=%s'%(kapdir,Y,Y15)
    cmd2 = '%s/cmult in=%s scalar=15.04 out=%s'%(kapdir,V,V17)
    os.system(cmd1)
    os.system(cmd2)
    cmd = '%s/add in1=%s in2=%s out=%s'%(kapdir,Y15,V17,YV)
    os.system(cmd)
    cmd = '%s/csub in=%s scalar=31.97 out=%s'%(kapdir,YV,dinominator)
    os.system(cmd)
    print '-Calculating the final uncertainty'
    cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,numerator,dinominator,fraction)
    os.system(cmd)
    cmd = '%s/mult in1=%s in2=%s out=%s'%(kapdir,fraction,ferror,Temp_error)
    os.system(cmd)
    cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,Temp_error,tempRAW,Frac_err)
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=100 out=%s'%(kapdir,Frac_err,Temp_percent)
    os.system(cmd)
    print '-Uncertainty calculation complete'

#################################################################################
########### Fixed variables
#Set FWHM for SCUBA-2 maps - MB (Main Beam) S (Secondary Beam)
fwhm450MB = 7.9 
fwhm850MB = 13.0 
fwhm450SB = 25.0 
fwhm850SB = 48.0 

sig450MB = FWHM(fwhm450MB,'TRUE')
sig850MB = FWHM(fwhm850MB,'TRUE')
sig450SB = FWHM(fwhm450SB,'TRUE')
sig850SB = FWHM(fwhm850SB,'TRUE')

#set normalisation constants
alpha450 = 0.94
alpha850 = 0.98
beta450 = 0.06
beta850 = 0.02

#Detection limit 
snr = 5

#################################################################################
########### free variables
#User selects a map input
#map = str(sys.argv[1])

#Hard Code directories you are working in !here!
#input_dir = 'input'
output_dir = 'output'
input_dir = str(sys.argv[4])

#enter the names of the SCUBA-2 data
file = str(sys.argv[1])
file450 = str(sys.argv[2])
file850 = str(sys.argv[3])

###################################################################
#Start of bulk script
###################################################################
############### 1 -- directory set up  #############
#cleans up old maps to be overwritten
print 1
#set up directories. 
if not os.path.exists(output_dir+'/map'):
    os.makedirs(output_dir+'/map')
if not os.path.exists(output_dir+'/map/'+file):
    os.makedirs(output_dir+'/map/'+file)
if not os.path.exists(output_dir+'/process'):
    os.makedirs(output_dir+'/process')
if not os.path.exists(output_dir+'/process/'+file):
    os.makedirs(output_dir+'/process/'+file)
if not os.path.exists(output_dir+'/math'):
    os.makedirs(output_dir+'/math')

######### 2 -- Convolve maps with the beam size of the other #########
print 2
#Convolution of maps
#4"/pixel for 450 micron pre-alignment; 6"/pixel for 850micron, convolve with other telescope beam

#get pixel size from image headers.
p450 = PARGET(input_dir+'/'+file450+'.sdf','fpixscale','ndftrace')
p850 = PARGET(input_dir+'/'+file850+'.sdf','fpixscale','ndftrace')
print "Pixel sizes are "+str(p450[0])+" and "+str(p850[0])
p450 = p450[0]
p850 = p850[0]
#Primary beam convolution using FWHM taken from Demspey2013.

print "primary beam convolution"
smooth450MB = fwhm850MB/p450
smooth850MB = fwhm450MB/p850
cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, input_dir+'/'+file450+'.sdf', smooth450MB, output_dir+'/process/'+file+'/s450convolveMB.sdf')
os.system(cmd)
cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, input_dir+'/'+file850+'.sdf', smooth850MB, output_dir+'/process/'+file+'/s850convolveMB.sdf')
os.system(cmd)

print "seconary beam convolution"
smooth450SB = fwhm850SB/p450
smooth850SB = fwhm450SB/p850
cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, input_dir+'/'+file450+'.sdf', smooth450SB, output_dir+'/process/'+file+'/s450convolveSB.sdf')
os.system(cmd)
cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, input_dir+'/'+file850+'.sdf', smooth850SB, output_dir+'/process/'+file+'/s850convolveSB.sdf')
os.system(cmd)

#Deffine sigma per pix size for use in normalisation
SIG_84_MB = sig850MB/p450
SIG_48_MB = sig450MB/p850
SIG_84_SB = sig850SB/p450
SIG_48_SB = sig450SB/p850

 #A (and B) = 1/(2*pi*(a4*(FM4**2)+b4*(FS4**2))) or 1/(2*pi*(a8*(FM8**2)+b8*(FS8**2)))
A = 1/(2.*np.pi*(alpha450*(SIG_48_MB**2.)+beta450*(SIG_48_SB**2.))) 
B = 1/(2.*np.pi*(alpha850*(SIG_84_MB**2.)+beta850*(SIG_84_SB**2.))) 

#set numbers for multiplication with convovle maps - in general they are of the form 2pi*(sig/pix)^2*A*alpha
A_a = 2.*np.pi*alpha450*A*(SIG_48_MB**2.)
A_b = 2.*np.pi*beta450*A*(SIG_48_SB**2.)
B_a = 2.*np.pi*alpha850*B*(SIG_84_MB**2.)
B_b = 2.*np.pi*beta850*B*(SIG_84_SB**2.)

if (float('%.f'%(A_a+A_b)) == float('%.f'%(B_a+B_b))) == 1.0:
    print 'Normalisation OK'
    print 'A = ',A
    print 'B = ',B
elif(float('%.f'%(A)) <> float('%.f'%(B))):
    print 'Normalistion FAIL'
    print 'A = ',A
    print 'B = ',B

print 'printing normalisation constants'
print '450: ',A_a,'+',A_b,' = ',A_a+A_b
print '850: ',B_a,'+',B_b,' = ',B_a+B_b

cmd1 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/process/'+file+'/s450convolveMB.sdf',B_a, output_dir+'/process/'+file+'/s450normMB.sdf')
cmd2 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/process/'+file+'/s450convolveSB.sdf',B_b, output_dir+'/process/'+file+'/s450normSB.sdf')
cmd3 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/process/'+file+'/s850convolveMB.sdf',A_a, output_dir+'/process/'+file+'/s850normMB.sdf')
cmd4 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/process/'+file+'/s850convolveSB.sdf',A_b, output_dir+'/process/'+file+'/s850normSB.sdf')
os.system(cmd1)
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)
cmd1 = '%s/add in1=%s in2=%s out=%s'%(kapdir,output_dir+'/process/'+file+'/s450normMB.sdf',output_dir+'/process/'+file+'/s450normSB.sdf', output_dir+'/process/'+file+'/s450convolve.sdf')
cmd2 = '%s/add in1=%s in2=%s out=%s'%(kapdir,output_dir+'/process/'+file+'/s850normMB.sdf',output_dir+'/process/'+file+'/s850normSB.sdf', output_dir+'/process/'+file+'/s850convolve.sdf')
os.system(cmd1)
os.system(cmd2)

############### 3 -- Align the 450 map with the 850 map ############
print 3
#Align 450 onto 850 convolved maps (in 2 or 3D where necessary)


percent = 5

SNR = 5

STDV_450 = noise.noise_by_data(input_dir+'/'+file450+'.sdf','FLASE')
STDV_850 = noise.noise_by_data(input_dir+'/'+file850+'.sdf','FLASE')

if STDV_850 == 0.0:
    TDV_850 = STDV_450/5.5
    print 'No noise for 850um. Value estimated from ratio and 450 level'

sigma450 = STDV_450*SNR
sigma850 = STDV_850*SNR

var450 = (STDV_450**2.0)*((alpha450**2.0) + (beta450**2.0))*(2./3.)
var850 = (STDV_850**2.0)*((alpha850**2.0) + (beta850**2.0))*(2./3.)


#Creating mask of input maps at 5sigma - for masking Align/collapse maps with
cmd1 = '%s/thresh in=%s out=%s thrlo=%f newlo=%s thrhi=%f newhi=%f QUIET'%(kapdir,input_dir+'/'+file450+'.sdf',output_dir+'/process/'+file+'/s450mask.sdf',sigma450,'bad',sigma450,1)
cmd2 = '%s/thresh in=%s out=%s thrlo=%f newlo=%s thrhi=%f newhi=%f QUIET'%(kapdir,input_dir+'/'+file850+'.sdf',output_dir+'/process/'+file+'/s850mask.sdf',sigma850,'bad',sigma850,1.)
os.system(cmd1)
os.system(cmd2)

print "Alignment of maps"
#get number of dimensions for maps to check if itermap cut or not
ndim450 = PARGET(output_dir+'/process/'+file+'/s450convolve.sdf','ndim','ndftrace')
ndim850 = PARGET(output_dir+'/process/'+file+'/s850convolve.sdf','ndim','ndftrace')
print 'Dimensions are: 450 (',ndim450,') + 850 (',ndim850,')'

#Sort out get parameters later.
if ndim450 == ndim850 == 3:
    print "path: 3 3"
    cmd1 = '%s/collapse in=%s axis=3 out=%s QUIET'%(kapdir,output_dir+'/process/'+file+'/s850convolve.sdf',output_dir+'/process/'+file+'/s850collapse.sdf')
    cmd2 = '%s/collapse in=%s axis=3 out=%s QUIET'%(kapdir,output_dir+'/process/'+file+'/s450convolve.sdf',output_dir+'/process/'+file+'/s450collapse.sdf')
    cmd3 = '%s/collapse in=%s axis=3 out=%s QUIET'%(kapdir,output_dir+'/process/'+file+'/s450mask.sdf',output_dir+'/process/'+file+'/s450maskcollapse.sdf')
    os.system(cmd1)
    os.system(cmd2) 
    os.system(cmd3)
    cmd1 = '%s/wcsalign in=%s out=%s ref=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,output_dir+'/process/'+file+'/s450collapse.sdf',output_dir+'/process/'+file+'/s450align.sdf',output_dir+'/process/'+file+'/s850collapse.sdf')
    cmd2 = '%s/wcsalign in=%s out=%s ref=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,output_dir+'/process/'+file+'/s450maskcollapse.sdf',output_dir+'/process/'+file+'/s450maskalign.sdf',output_dir+'/process/'+file+'/s850collapse.sdf')
    os.system(cmd1)
    os.system(cmd2)
else:
    print 'maps are not in 3D - script will Break here'
    

############## 4-Thresh maps  ############
print 4
#Threshes maps to remove negative components from scaled maps, also clips maps based on SNR 
#first - divide the 450 mask by 2.25 to account for align bias.
align = (p850**2.)/(p450**2.)
print 'Align = '+str(align)

print '5sigma noise level:'
print '450 = '+str(sigma450)
print '850 = '+str(sigma850)

cmd = '%s/cdiv in=%s scalar=%f out=%s'%(kapdir,output_dir+'/process/'+file+'/s450maskalign.sdf',align,output_dir+'/process/'+file+'/s450readymask.sdf')
os.system(cmd)
cmd1 = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,output_dir+'/process/'+file+'/s450align.sdf',output_dir+'/process/'+file+'/s450readymask.sdf',output_dir+'/process/'+file+'/s450alignTH.sdf')
cmd2 = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,output_dir+'/process/'+file+'/s850collapse.sdf',output_dir+'/process/'+file+'/s850mask.sdf',output_dir+'/process/'+file+'/s850collapseTH.sdf')
os.system(cmd1)
os.system(cmd2)

############### 5 -- Make Ratio maps ###########################
print 5
#checks for variance array for maps to allow snr cuts
### IF statement on when to make cuts (no. of dimensions and wether there is varience) ###
print 'Calculating Flux Ratio, stored in: '+output_dir+'/map/'+file+'/'+file+'_Sratio.sdf'
print 'Number of dimensions: ',ndim850
cmd = '%s/div in1=%s in2=%s out=%s'%(kapdir,output_dir+'/process/'+file+'/s450alignTH.sdf',output_dir+'/process/'+file+'/s850collapseTH.sdf',output_dir+'/map/'+file+'/'+file+'_Sratio.sdf')
os.system(cmd)

Sratio = str(output_dir+'/map/'+file+'/'+file+'_Sratio.sdf')

######### 6 -- creates new uncertainty arrays #########
print 6
# creates new Varience array as existing one is not correct (exact reason unknown)
print 'Creating new uncertainty array of ratio  maps:'
#Error calculated by adding fractonal errors in quadrature.

n450 = output_dir+'/process/'+file+'/s450noise.sdf' #maps of constant noise level
n850 = output_dir+'/process/'+file+'/s850noise.sdf' #maps of constant noise level
s450 = output_dir+'/process/'+file+'/s450align.sdf' #maps of constant flux
s850 = output_dir+'/process/'+file+'/s850collapse.sdf' #maps of constant flux
f450 =  output_dir+'/process/'+file+'/s450frac.sdf' #maps of fractional error
f450sq =  output_dir+'/process/'+file+'/s450sqfrac.sdf' #maps of fractional error Squared
f850 =  output_dir+'/process/'+file+'/s850frac.sdf' #maps of fractional error
f850sq =  output_dir+'/process/'+file+'/s850sqfrac.sdf' #maps of fractional error Squared
fsum = output_dir+'/process/'+file+'/sumfrac.sdf' #maps of sum of frac errors
ferror = output_dir+'/process/'+file+'/fracerror.sdf' #maps of frac error

error = output_dir+'/map/'+file+'/'+file+'Sratio_error.sdf' #maps of error in ratio
var = output_dir+'/map/'+file+'/'+file+'Sratio_var.sdf' #maps of error in ratio

cmd1 = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,s450,n450,0,0,var450,var450)
cmd2 = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%f newhi=%f QUIET'%(kapdir,s850,n850,0,0,var850,var850)
os.system(cmd1)
os.system(cmd2)
cmd1 = '%s/div in1=%s in2=%s out=%s'%(kapdir,n450,s450,f450)
cmd2 = '%s/div in1=%s in2=%s out=%s'%(kapdir,n850,s850,f850)
os.system(cmd1)
os.system(cmd2)
cmd1 = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,f450,f450,f450sq)
cmd2 = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,f850,f850,f850sq)
os.system(cmd1)
os.system(cmd2)
cmd = '%s/add in1=%s in2=%s out=%s'%(kapdir,f450sq,f850sq,fsum)
os.system(cmd)
cmd = '%s/pow in=%s power=0.5 out=%s'%(kapdir,fsum,ferror)
os.system(cmd)
cmd = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,ferror,Sratio,error)
os.system(cmd)
cmd = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,error,error,var)
os.system(cmd)

######### 7 -- Create FITS files and pass to 'temperature.py' ###########
print 7
#Produces RAW temperature map with no tidying or corrections due to variance.
type = 't' # str(sys.argv[5])
while True:
    if str(type) == 't':
    #creates FITS file to output to python
        ratio_fit = str(output_dir+'/process/'+file+'/'+file+'Sratio.fit')
        var_fit = str(output_dir+'/process/'+file+'/'+file+'Sratio_var.fit')
        tempRAW = str(output_dir+'/process/'+file+'/'+file+'_tempRAW.sdf')
        
        cmd = 'rm %s'%(ratio_fit)
        os.system(cmd)
        cmd = 'rm %s'%(var_fit)
        os.system(cmd)
        cmd1 = '%s/ndf2fits in=%s out=%s QUIET'%(convdir,Sratio,ratio_fit)
        cmd2 = '%s/ndf2fits in=%s out=%s QUIET'%(convdir,var,var_fit)
        os.system(cmd1)
        os.system(cmd2)

        beta = 1.8 #float(sys.argv[6])
    #Move to Temperature.py
        temperature(Sratio,ratio_fit,var_fit,beta,'FALSE')

        cmd1 = '%s/fits2ndf in=%s out=%s QUIET'%(convdir,ratio_fit,tempRAW)
        os.system(cmd1)
        break
    elif str(type) == 'b':
        T = 15. #float(sys.argv[6])
        betamap = output_dir+'/map/'+file+'/'+file+'_betaRaw.sdf'
        beta = 'log(IA*(exp(31.97/'+str(T)+')-1)/(exp(16.93/'+str(T)+')-1))/log(450/850))-3'
        cmd = '%s/maths %s ia=%s out=%s'%(kapdir,beta,Sratio,betamap)
        os.system(cmd)
        break
    type = str(raw_input("Please enter 'b' or 't':")) #prompt user to enter correct value if they get it wrong

######### 8 -- Makes Variance arrays for Temp #########
print 8
# creates new Varience array for temp. maps from Analytical solution 
print 'Creating new error array of temp  maps:'
tempCUT = str(output_dir+'/process/'+file+'/'+file+'_tempcut.sdf')
Temp_error =  output_dir+'/map/'+file+'/'+file+'temp_error.sdf' #maps of error in temp
Temp_percent =  output_dir+'/process/'+file+'/'+file+'temp_errorPC.sdf' #maps of percentage error in temp.
Temp_mask = output_dir+'/process/'+file+'/'+file+'temp_'+str(percent)+"mask.sdf" #maps of percentage error in temp.
Temp_mask1 =  output_dir+'/math/'+file+'temp_'+str(percent)+"mask1.sdf" #maps of percentage error in temp.
Temp_mask2 =  output_dir+'/math/'+file+'temp_'+str(percent)+"mask2.sdf" #maps of percentage error in temp.
TempFINAL = output_dir+'/map/'+file+'/'+file+"temperature.sdf"
#cuts anything calculated > 1000K, as unphysical.
print "Cutting high (>999.98K) pixels:"
cmd = '%s/thresh in=%s out=%s thrlo=0 newlo=0 thrhi=%f newhi=%s QUIET'%(kapdir,tempRAW,tempCUT,999.98,'bad')
os.system(cmd)

#Function to calculate maps of temperature uncertainity
print 'Calculating error in temperature'
temp_error(Sratio,tempCUT,ferror,Temp_error,Temp_percent)

print 'Precent cut at: ',percent,'%'
#Cut anything with an uncertainty higher than a given percent
cmd = '%s/thresh in=%s out=%s thrlo=%f newlo=%f thrhi=%f newhi=%s QUIET'%(kapdir,Temp_percent,Temp_mask1,percent,1,percent,'bad') #cutting where error > 5%
os.system(cmd)
cmd = '%s/thresh in=%s out=%s thrlo=%f newlo=%s thrhi=%f newhi=%s QUIET'%(kapdir,Temp_percent,Temp_mask2,0,'bad',0,1) #cutting where error is less than 0 (don't know why this is occouring)
os.system(cmd)
cmd = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,Temp_mask1,Temp_mask2,Temp_mask)
os.system(cmd)
cmd = '%s/mult in1=%s in2=%s out=%s QUIET'%(kapdir,tempCUT,Temp_mask,TempFINAL)
os.system(cmd)

print 'Create temperature map: ',TempFINAL

print '========================'
cmd = '%s/stats ndf=%s'%(kapdir,TempFINAL)
os.system(cmd)

cmd1 = 'rm -r %s/math/'%(output_dir)
os.system(cmd1)
print 'End'
