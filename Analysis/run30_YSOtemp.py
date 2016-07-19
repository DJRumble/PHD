import numpy as np
from pylab import *
from scipy import stats
import matplotlib.pyplot as plt
import time

data =  np.loadtxt('data/YSOtemp_match.txt',dtype='float',comments="#")

Ptemp0 = []
dPtemp0 = []
Ctemp0 = []
dCtemp0 = []
Ptemp1 = []
dPtemp1 = []
Ctemp1 = []
dCtemp1 = []
PtempFS = []
dPtempFS = []
CtempFS = []
dCtempFS = []

#protostar temp
for i in data:
    if i[7] <= 70.:
        Ptemp0.append(i[2])
        dPtemp0.append(i[3])
        Ctemp0.append(i[4])
        dCtemp0.append(i[5])
    elif 350. >= i[7] >= 70.:
        Ptemp1.append(i[2])
        dPtemp1.append(i[3])
        Ctemp1.append(i[4])
        dCtemp1.append(i[5])
    elif i[7] >= 350.:
        PtempFS.append(i[2])
        dPtempFS.append(i[3])
        CtempFS.append(i[4])
        dCtempFS.append(i[5])
    else:
        print 'No Tbol info'

##########
###PLOT###
##########
xlabel('Envelope temp. (K)')
ylabel('Clump temp. (K)')
#C0
errorbar(Ptemp0,Ctemp0,xerr=dPtemp0,yerr=dCtemp0,fmt='go',label='Class 0')
#C1
errorbar(Ptemp1,Ctemp1,xerr=dPtemp1,yerr=dCtemp1,fmt='bo',label='Class I')
#CII
errorbar(PtempFS,CtempFS,xerr=dPtempFS,yerr=dCtempFS,fmt='ro',label='Class FS')

ylim([5,39])
xlim([5,39])
#plt.axvline(x=1.2,linestyle='dashed',color='red')

legend(loc=4)

A = [0,40]

plot(A,A,'k-')

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_YSOtemp.pdf'%(date))

show()
