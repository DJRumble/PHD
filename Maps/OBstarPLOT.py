import numpy as np
from pylab import *
import matplotlib.pyplot as plt

OB = np.loadtxt('data/obcatalog_deg.txt',dtype='float')
YSO = np.loadtxt('data/GBS_YSO_master_FLT.tab',dtype='float')
Clumps = np.loadtxt('SMM/SMM_master_FLT.tab',dtype='float')
OB_GBS = np.loadtxt('OBstars_list.txt',dtype='float')

x_OB = OB[:,1]
y_OB = OB[:,2]

x_YSO = YSO[:,0]
y_YSO = YSO[:,1]

x_Clumps = Clumps[:,0]
y_Clumps = Clumps[:,1]

x_GBS = OB_GBS[:,0]
y_GBS = OB_GBS[:,1]


plt.scatter(x_OB,y_OB,marker='o',color='r',s=5)
plt.scatter(x_Clumps,y_Clumps,color='b',s=6)
plt.scatter(x_YSO,y_YSO,marker='+',facecolor='k',s=200)
plt.scatter(x_GBS,y_GBS,marker='*',edgecolors='k',color='w',s=200)


plt.annotate('Cepheus', xy=(320, 82),color='k')
plt.annotate('IC5146', xy=(320, 50),color='k')
plt.annotate('Ophiuchus', xy=(240, -25),color='k')
plt.annotate('Serpens-Aquila', xy=(265, 0),color='k')
plt.annotate('Lupus', xy=(230, -40),color='k')
plt.annotate('Corona Australis', xy=(320, -45),color='k')
plt.annotate('Musca', xy=(180, -70),color='k')
plt.annotate('Chamaeleon', xy=(175, -85),color='k')
plt.annotate('Orion', xy=(80, -12),color='k')
plt.annotate('Perseus', xy=(45, 30),color='k')
plt.annotate('Auriga', xy=(90, 40),color='k')
plt.annotate('Taurus', xy=(60, 20),color='k')


plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')

xlim([360,0])
ylim([-90,90])



plt.show()
