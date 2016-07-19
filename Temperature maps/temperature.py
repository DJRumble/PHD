import numpy as np
import pyfits
import sys

####### Function A ##################
def zoominratio(min, max, count):
    mid = (min+max)/2  #sets middle
    count += 1         #increments flag
    #print count
    #print mid
    if(count < 14):
        #print lookup[mid][1]
        #print ratio[i][j]
        #if value in top half of search, search again with top half of section of array
        if(ratio[i][j] > lookup[mid][1]):
            zoominratio(int(mid)-1, max, count)
        #if in bottom half, call bottom half of array again
        elif(ratio[i][j] < lookup[mid][1]): 
            zoominratio(min, int(mid)+1, count)
    #once zoomed into small enough section of array, find best fit
    else: 
        result = 1                                                             
        difference = 100
        for k in range(min, max):
            #check if closest match or not
            if (np.abs(lookup[k][1] - ratio[i][j]) < difference):
                result = k #set array index if smallest so far
                difference = np.abs(lookup[k][1] - ratio[i][j]) #update smallest difference 
            #print k
            #print result
        ratio[i][j] = lookup[result][0] 
    #output correct temperature once best fit found
    return

###### Function B ###################
#same as zoominratio, just for upper limit
def zoominupper(min, max, count):
    mid = (min+max)/2
    count += 1
    if(count < 14):
        if(upperratio[i][j] > lookup[mid][1]):
            zoominupper(int(mid)-1, max, count)
        elif(upperratio[i][j] < lookup[mid][1]):
            zoominupper(min, int(mid)+1, count)
    else:    
        result = 1
        difference = 100
        for k in range(min, max):
            if (np.abs(lookup[k][1] - upperratio[i][j]) < difference):          
                result = k
                difference = np.abs(lookup[k][1] - upperratio[i][j])
        upperratio[i][j] = lookup[result][0]
    return

###### Function C ################
#same as zoominratio, just lower limit
def zoominlower(min, max, count):
    mid = (min+max)/2
    count += 1
    if(count < 14):
        if(lowerratio[i][j] > lookup[mid][1]):
            zoominlower(int(mid)-1, max, count)
        elif(lowerratio[i][j] < lookup[mid][1]):
            zoominlower(min, int(mid)+1, count)
    else:        
        result = 1
        difference = 100
        for k in range(min, max):                                               
            if (np.abs(lookup[k][1] - lowerratio[i][j]) < difference):          
                result = k
                difference = np.abs(lookup[k][1] - lowerratio[i][j])
        lowerratio[i][j] = lookup[result][0]               
    return


##### DEFFINE Parrameters ##########
beta = float(sys.argv[5])     #beta value held constant for calculations
row = int(sys.argv[4])        #dimensions of table
column = int(sys.argv[3])
isvar = sys.argv[6]           #variance flag
image = pyfits.open('{0}/{1}.fit'.format(sys.argv[1], sys.argv[2]), mode = 'update')
ratio = image[0].data         #opens FITS file created from NDF
if(isvar == "yesvar"):
    variance = image[1].data  #opens variance array if NDF has one
lookup = np.zeros((99501, 3), float)
upperratio = np.zeros((row, column), float)
lowerratio = np.zeros((row, column), float)
T = 5.0

##### Main Program #################

print 'Beta = ' + str(beta)
print 'row = ' + str(row)
print 'column =' + str(column)
print 'isvar =' + str(isvar)
print 'image =' + str(image)

#For range of temp. of  from 5.0 going up in increments of 0.01, specific ratios are calculated.
for i in range(0, 99501):
    lookup[i][0] = float(T)
    T += 0.01 #creates temperature-ratio lookup table
    lookup[i][1] = (17**(3+beta)) / (9**(3+beta)) * (np.exp(16.93/lookup[i][0])-1) / (np.exp(31.97/lookup[i][0])-1)

#compare table of ratios with real data and use this to build new temp maps - based on varience
for i in range(0, row):
    for j in range(0, column):
        #if variance flag set, find upper and lower limits based on ratio variance
        if (ratio[i][j] > 0 and isvar == "yesvar"):
            upperratio[i][j] = ratio[i][j] + variance[i][j]
            lowerratio[i][j] = ratio[i][j] - variance[i][j]

for i in range(0, row):
    for j in range(0, column):
        if ratio[i][j] > 0:
            zoominratio(0, 99500, 0) #call functions that find best temperature fit for ratios
           # print ratio[i][j] 
            if(isvar == "yesvar"):
                zoominupper(0, 99500, 0)
                zoominlower(0, 99500, 0) #call limits if applicable and calculate variance by difference of limits
                variance[i][j] = (upperratio[i][j] - lowerratio[i][j])/2

#output new values (temperature from ratio) back to FITS file and save
image.flush() 
image.close()           
       
           
