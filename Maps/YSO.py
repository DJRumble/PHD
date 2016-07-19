"""Module for plot YSO from the various sources - used in plotting figures"""

import atpy
import c2dtable
import numpy 
import astropy.io.fits as pyfits

#key:
#Protostars - Green
#Class 0  
#Class I 
#Class FS 

#Pre-Main-Sequence Stars - red
#Class II 
#Cass TD
#Class III 

def GBS(fig):
    #Plots from the GBS Spitzer catalog
    #Input a figure in '.sdf' format and this module will arrange for some markers to plotted over thier location. To be used in conjunction with the Map Maker

####################################################
#load GBS catalogue 

    tbl = atpy.Table()
    c2dtable.c2d_read('catalog-Serp-YSOc+head.tbl',tbl,3)
    ra_t=tbl.data['RA']
    dec_t=tbl.data['DEC']
    alpha= tbl.data['alpha']

###################################################
#Classification by Colour

    CYSO = []
    SYSO = []
    FYSO = []
    i = 0

    while (i < len(alpha)): #iterate through length of the array   
        a = alpha[i]
        if a >= 0.3: #Class I as yellow
            colour = 'lime'
            shape = '^'
            fill = 'lime'
        elif 0.3 > a >= -0.3: #Class FS as Orange
            colour = 'lime'
            shape = 's'
            fill = 'lime'
        elif -0.3 > a >= -1.6: #Class II as Red
            colour = 'red'
            shape = 'o'
            fill = colour
        elif -1.6 > a: #Class III as Magenta
            colour = 'red'
            shape = 'o'
            fill = 'none'
        else: #else(Clss O)
            colour = 'lime'
            shape = 'x'
            fill = colour
        CYSO.append(colour)
        SYSO.append(shape)
        FYSO.append(fill)
        i = i + 1
    
    ra = 200.0+ra_t
    dec = dec_t
    fig.show_markers(ra,dec, layer='GBS',edgecolor=CYSO ,facecolor=FYSO, marker='o',s=40, alpha=1)   
    return

def c2dGBS(fig):
    #Plots from the c2d+GBS YSO catalogue from Dunham 2015 + additions
    #Input a figure in '.sdf' format and this module will arrange for some markers to plotted over thier location. To be used in conjunction with the Map Maker

####################################################
#load c2d+GBS catalogue 

    data = numpy.loadtxt('/data/damian/run27/data/GBS_YSO_master.tab',dtype='string')
    ra = (numpy.array(data[:,2])).astype(numpy.float)
    dec = (numpy.array(data[:,3])).astype(numpy.float)
    alpha = (numpy.array(data[:,4])).astype(numpy.float)
    tbol = (numpy.array(data[:,5])).astype(numpy.float)

###################################################
#Classification by Colour

    CYSO = []
    SYSO = []
    FYSO = []
    i = 0

    while (i < len(alpha)): #iterate through length of the array   
        a = alpha[i]
        b = tbol[i]
        if 0.3 > a >= -0.3: #Class FS
            colour = 'red'
            fill = 'lime'
        elif -0.3 > a >= -1.6: #Class II 
            colour = 'k'
            fill = 'red'
        elif -1.6 > a: #Class III
            colour = 'red'
            fill = 'none'
        if a >= 0.3: 
            if b <= 70.: #Class 0 
                colour = 'lime'
                fill = 'none'
            else: #Class I 
                colour = 'k' 
                fill = 'lime'

        CYSO.append(colour)
        FYSO.append(fill)
        i = i + 1

    fig.show_markers(ra,dec, layer='GBS',edgecolor=CYSO ,facecolor=FYSO, marker='o',s=25, alpha=1)   
    #fig.show_markers(ra,dec, layer='GBS2',edgecolor=CYSO ,facecolor=FYSO, marker='x',s=120, alpha=1)  
    return

def c2dGBS_v2(fig):
    #Plots from the c2d+GBS YSO catalogue from Dunham 2015 + additions
    #Input a figure in '.sdf' format and this module will arrange for some markers to plotted over thier location. To be used in conjunction with the Map Maker

####################################################
#load c2d+GBS catalogue 

    data = numpy.loadtxt('/data/damian/run27/data/GBS_YSO_master.tab',dtype='string')

    ra = (numpy.array(data[:,2])).astype(numpy.float)
    dec = (numpy.array(data[:,3])).astype(numpy.float)
    alpha = (numpy.array(data[:,4])).astype(numpy.float)
    tbol = (numpy.array(data[:,5])).astype(numpy.float)

###################################################
#Classification by Colour

    raC0,decC0 = [],[]
    raCI,decCI = [],[]
    raCII,decCII = [],[]
    raCIII,decCIII = [],[]

    i = 0
    while (i < len(alpha)): #iterate through length of the array   
        a = alpha[i]
        b = tbol[i]
        if 0.3 > a >= -0.3: #Class FS
            raCI.append(ra[i])
            decCI.append(dec[i])
        elif -0.3 > a >= -1.6: #Class II 
            raCII.append(ra[i])
            decCII.append(dec[i])
        elif -1.6 > a: #Class III
            raCIII.append(ra[i])
            decCIII.append(dec[i])
        if a >= 0.3: 
            if b <= 70.: #Class 0 
                raC0.append(ra[i])
                decC0.append(dec[i])
            else: #Class I 
                raCI.append(ra[i])
                decCI.append(dec[i])
        i = i + 1

    fig.show_markers(raC0,decC0,edgecolor='g' ,facecolor='none', marker=('o'),s=250, alpha=1)  
    fig.show_markers(raCI,decCI,edgecolor='g' ,facecolor='g', marker=('+'),s=250, alpha=1) 
    fig.show_markers(raCII,decCII,edgecolor='r' ,facecolor='r', marker=('x'),s=250, alpha=1)  
    fig.show_markers(raCIII,decCIII,edgecolor='r' ,facecolor='none', marker=('*'),s=250, alpha=1) 
    return


def damiani(fig):
     #plots position of X-ray emission partially covering MWC297

####################################################
#Load mwc297 data catalogs

     a = numpy.loadtxt('/data/damian/data/20130125_d06_mwc297_xraysources_deg.txt',dtype='float')
     ra = a[:,0]
     dec = a[:,1]
     fig.show_markers(ra,dec, layer='damiani',edgecolor='blue' ,facecolor='none', marker='x',s=40, alpha=1)
     return

def w40_radio(fig):
    #Plots OBstars from literature unique to W40. 

####################################################
#Load W40 data catalog

    #data = numpy.loadtxt('/data/damian/data/aquila/UCHIIc_FW_W40.tab')
    #ra = data[:,1]
    #dec = data[:,2]

    data = numpy.loadtxt('/data/damian/data/aquila/Ortiz-Leon-radio-W40-YSO_flt.tab')
    ra = data[:,7]
    dec = data[:,8]
    fig.show_markers(ra,dec, layer='radio2',edgecolor='blue' ,facecolor='none', marker='o',s=70, alpha=1)

    data = numpy.loadtxt('/data/damian/data/aquila/Ortiz-Leon-Rodriguez_YSOmatches.tab')
    ra = data[:,0]
    dec = data[:,1]
    fig.show_markers(ra,dec, layer='radio4',edgecolor='blue' ,facecolor='blue', marker='o',s=70, alpha=1)


    ra = [277.84071,277.83919,277.84070,277.84071]
    dec = [-2.106645,-2.104122,-2.112197,-2.116741]
    fig.show_markers(ra,dec, layer='radio3',edgecolor='blue' ,facecolor='blue', marker='x',s=70, alpha=1)

    #w40 = numpy.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    #ra = w40[:,0]
    #dec = w40[:,1]
    #fig.show_markers(ra,dec, layer='OB obs',edgecolor='black' ,facecolor='yellow', marker='o',s=40, alpha=1)

def w40(fig):
    #Plots OBstars from literature unique to W40. 

####################################################
#Load W40 data catalogs

    ysos = numpy.loadtxt('/data/damian/data/aquila/W40_YSOs_biglist.tab')
    ra = ysos[:,0]
    dec = ysos[:,1]
    obj = ysos[:,2]

    CYSO = []
    SYSO = []
    FYSO = []
    i = 0

    while (i < len(ysos)): #iterate through length of the array   
        a = obj[i]
        if a == 1: #Class I 
            colour = 'k' 
            fill = 'lime'
        elif a == 2: #Class II 
            colour = 'k'
            fill = 'red'
        elif a == 3: #Class III 
            colour = 'red'
            fill = 'none'
        elif a == 4: #Class III 
            colour = 'black'
            fill = 'yellow'
        elif a == 0: #Class 0
            colour = 'lime'
            fill = 'none'
        CYSO.append(colour)
        FYSO.append(fill)
        i = i + 1
    fig.show_markers(ra,dec, layer='biglist',edgecolor=CYSO ,facecolor=FYSO, marker='o',s=40, alpha=1)  

    #w40 = numpy.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    
    #ra = w40[:,0]
    #dec = w40[:,1]

    #fig.show_markers(ra,dec, layer='OB obs',edgecolor='black' ,facecolor='yellow', marker='o',s=40, alpha=0.5)

    return


def mwc297(fig):
    #Plots sources from literature unique to MWC 297. No specification for source all files are within one file.

####################################################
#Load mwc297 data catalogs
    CYSO = []
    SYSO = []
    FYSO = []
    a = numpy.loadtxt('/data/damian/data/20140219_MWC297_YSO.txt',dtype='float')
    ra = a[:,0]
    dec = a[:,1]
    alpha = a[:,2]
    i = 0
    while (i < a.shape[0]): #iterate through length of the array
        if alpha[i] >= 0.3: #Class I as yellow
            facecolour = 'lime'
            shape = '^'
            fill = 'none'    
        elif 0.3 > alpha[i] >= -0.3: #Class FS as Orange
            colour = 'lime'
            shape = 's'
            fill =  'none'
        elif -0.3 > alpha[i] >= -1.6: #Class II as Red
            colour = 'red'
            shape = 'o'
            fill = colour
        elif -1.6 > alpha[i]: #Class III as Magenta
            colour = 'red'
            shape = 'o'
            fill = 'none'
        else: #else(Clss O)
            colour = 'lime'
            shape = 'x'
            fill = colour
        CYSO.append(colour)
        SYSO.append(shape)
        FYSO.append(fill)
        i = i + 1

    fig.show_markers(ra,dec, layer='mwc297',edgecolor=CYSO ,facecolor=FYSO, marker='o',s=40, alpha=1)  
    return


def oth(fig):
     #Plots sources from literature for other regions. 

##################################################
#Load All YSO data

#key:
#Displays Protostars (classes 0,I,FS) as Lime 
#Displays PMS stars (classes II,D and III) as Red

#Markers - 

#Connelly +
#gutermuth v 
#gutermuth (GBS as above) d 
#evans s
#winston o
#Kuhn x
#damiani x
#IRO d
#h2o h

#Kunn Xray - x
    k1 = '/data/damian/data/aquila/20130419_k10_protostars.txt'#note these are columns 1,2
    k2 = '/data/damian/data/aquila/20130419_k10_PMS.txt'#note these are columns 1,2

    YSO = numpy.loadtxt(k1)
    ra = YSO[:,1]
    dec = YSO[:,2]
    fig.show_markers(ra,dec, layer='k1', edgecolor='lime',facecolor='lime',marker='d',s=50, alpha=0.75)

    YSO = numpy.loadtxt(k2)
    ra = YSO[:,1]
    dec = YSO[:,2]
    fig.show_markers(ra,dec, layer='k2', edgecolor='red',facecolor='red',marker='d',s=50, alpha=0.75)

#Spitzer Young Clusters - v
    g2 = '/data/damian/data/gutermuth09/20130411_G09_PMS_all.txt'
    g1 = '/data/damian/data/gutermuth09/20130411_G09_protostars_all.txt'

    #YSO = numpy.loadtxt(g1)
    ra = YSO[:,0]
    dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='g1', edgecolor='red',facecolor='red',marker='^',s=50, alpha=0.75)

    #YSO = numpy.loadtxt(g2)
    ra = YSO[:,0]
    dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='g2', edgecolor='red',facecolor='red',marker='^',s=50, alpha=0.75)

#Damiani06 X-ray - x
    d2 = '/data/damian/data/mwc297/20130125_d06_mwc297_xraysources_deg.txt'
    YSO = numpy.loadtxt(d2)
    ra = YSO[:,0]
    dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='d2', edgecolor='red',facecolor='red',marker='o',s=50, alpha=0.75)

#Evans YSOs - s
    e2 = '/data/damian/data/evans09/20130411_E09_PMS.txt'  #note these are columns 1,2
    e1 = '/data/damian/data/evans09/20130411_E09_protostars.txt' #note these are columns 1,2

    YSO = numpy.loadtxt(e1)
    ra = YSO[:,1]
    dec = YSO[:,2]
    fig.show_markers(ra,dec, layer='e1', edgecolor='k',facecolor='lime',marker='o',s=50, alpha=0.75)

    YSO = numpy.loadtxt(e2)
    ra = YSO[:,1]
    dec = YSO[:,2]
    fig.show_markers(ra,dec, layer='e2', edgecolor='k',facecolor='red',marker='o',s=50, alpha=0.75)

#Winston YSOs - o  
    w2 = '/data/damian/data/winston07/20130411_W07_PMS.txt'
    w1 = '/data/damian/data/winston07/20130411_W07_protostars.txt'

    YSO = numpy.loadtxt(w1)
    ra = YSO[:,0]
    dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='w1', edgecolor='lime',facecolor='lime',marker='o',s=50, alpha=0.75)

    YSO = numpy.loadtxt(w2)
    ra = YSO[:,0]
    dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='w2', edgecolor='red',facecolor='red',marker='o',s=50, alpha=0.75)

#Conelly Class Is - o
    c1 = '/data/damian/data/connelly10/20130211c10_c_0_I_south_deg.txt'

    YSO = numpy.loadtxt(c1)
    ra = YSO[:,0]
    dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='c1', edgecolor='lime',facecolor='lime',marker='o',s=50, alpha=0.75)

#other (3) - MWC297 specific

    i3 = '/data/damian/data/mwc297/20130208_IRO_serpens_E1.txt'
    s3 = '/data/damian/data/mwc297/20130208_combined_serpens_E1.txt'
    x3 = '/data/damian/data/stars.txt'

    #YSO = numpy.loadtxt(i3)
    #ra = YSO[:,0]
    #dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='i3', edgecolor='magenta',facecolor='none',marker='d',s=50, alpha=0.75)

    #YSO = numpy.loadtxt(s3)
    #ra = YSO[:,0]
    #dec = YSO[:,1]
    #fig.show_markers(ra,dec, layer='s3', edgecolor='cyan',facecolor='none',marker='h',s=50, alpha=0.75)
    return

def mwc297_oth(fig):
    #Show all Other objects unique to MWC 297 (such as companions)
    x3 = '/data/damian/data/stars.txt'

    YSO = numpy.loadtxt(x3)
    ra = YSO[:,0]
    dec = YSO[:,1]
    fig.show_markers(ra,dec, layer='xa', edgecolor='black', facecolor='yellow', marker='v',s=100, alpha=1)

    YSO = numpy.loadtxt(x3)
    ra = YSO[:,0]
    dec = YSO[:,1]
    fig.show_markers(ra,dec, layer='xb', edgecolor='black', facecolor='yellow', marker='^',s=100, alpha=1)

#Sh2-62
    ra = 276.916666667
    dec = -3.83166666667
    fig.show_markers(ra,dec, layer='sh2-62', edgecolor='lime', marker='h',s=100, alpha=1)

#Companion 1 (+2.4'' +14'')
    ra = 276.924708333
    dec = -3.835055555556
    fig.show_markers(ra,dec, layer='companion1', edgecolor='black',facecolor='yellow', marker='d',s=100, alpha=1)

#Companion 1 (+1'' -19'')
    ra = 276.918875 
    dec = -3.8258888888
    fig.show_markers(ra,dec, layer='companion2', edgecolor='black',facecolor='yellow', marker='d',s=100, alpha=1)
    return




if __name__ == "__main__":
    GBS(fig)
    mwc297(fig)
    oth(fig)
    mwc297_oth(fig)
    w40(fig)
    w40_radio(fig)
