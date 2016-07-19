import os
import numpy as np
import string

maps = np.loadtxt('maps_central_pix.txt',dtype='string')

list = open("maps4fluxextractor.txt","w")
list = open("maps4fluxextractor.txt","a")

#loop through maps
for i in maps:
    #split INPUT file
    s850 = string.split(i,',')[0]
    #deffine region name
    name = string.split(s850,'/')[1]
    print name
    file = string.split(name,'_850')[0]
    print file
    region = string.split(name,'_20')[0]
    print region

    temp = 'input_temp/'+str(file)+'-autotemperature.sdf'
    temp_err = 'input_temp/error/'+str(file)+'-autotemp_error.sdf'
    alpha = 'input_temp/alpha/'+str(file)+'-auto_Salpha.sdf'
    proto = 'YSO/points/'+str(region)+'_protostar.sdf'
    surface = 'YSO/surface/'+str(region)+'_YSOsurface_r300.sdf'
    clumps = 'FW/clumps/'+str(file)+'_clumps.sdf'
    cat = 'FW/cat/'+str(file)+'_clump_cat.fit'

    list.write('########################################################################################################################\n')
    list.write(s850+','+temp+','+temp_err+','+alpha+','+proto+','+surface+','+clumps+','+cat+'\n')

    #break
