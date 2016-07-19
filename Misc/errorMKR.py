#DJR UoE
#02/07/2013 - errorMKR.csh
#
#This code produces an error array for ~450850.sdf (ratio) maps
#
#there are three fundamental process:
#1) error in 450 nd 850 due to convolution (alpha and beta)
#2) error in alighnment, 450 only
#3) Final error in ratio. Division of align450 and collapse850

#Output files go to directory output_5
##################################################################

import pyfits

###################################################################
#Start of parameters
###################################################################




#STDV of maps - updated 17 june 2013 - http://wiki.astro.ex.ac.uk/bin/view/JCMTGouldBelt/SCUBA-2IR1NoiseEstimates
STDV_south_450 = 0.000194042513696 
STDV_south_850 = 0.0000475598756647
STDV_main_450 = 0.000328060247248 
STDV_main_850 = 0.0000776023912041 
STDV_mwc297_450 = 0.00176626163015  
STDV_mwc297_850 = 0.000220099947459 
STDV_E1_450 = 0.00202120458139 
STDV_E1_850 = 0.000187665752113 

###################################################################
#Start of functiuons
###################################################################

'''
#names of files
Aquila_20130528_s450_IR1_JH	
SerpensE1_20130425_s450_IR1_JH
SerpensMain_20130425_s450_IR1_JH
SerpensMWC297_20130425_s450_IR1_JH
'''

def convolve450(STDV_450):
    alpha450 = 0.94
    beta450 = 0.06
    
    convolve_450 = ((STDV_450**2.0)*((alpha450**2.0)+(beta450**2.0)))**(0.5)
    print convolve_450
    return 

def convolve850(STDV_850):
    alpha850 = 0.98
    beta850 =  0.02

    convolve_850 = ((STDV_850**2.0)*((alpha850**2.0)+(beta850**2.0)))**(0.5)
    print convolve850
    return

def align450(STDV_450):
    align450 = (2.0/3.0)*STDV_450
    print align450
    return 

s450align = pyfits.getdata("output_5/450/Aquila_20130528_s450_IR1_JH/s450align.sdf")
s850collapse = pyfits.getdata("output_5/850/Aquila_20130528_s450_IR1_JH/s850collapse.sdf")
ratio450850 = pyfits.getdata("output_5/map/Aquila_20130528_s450_IR1_JH/Aquila_20130528_IR1_JH450850.sdf")

A = (align450(convolve450(STDV_south_450))/s450align)**2.0
B = (convolve850(STDV_south_850)/s850collapse)**2.0

array = ratio450850*((A+B)**(0.5))


