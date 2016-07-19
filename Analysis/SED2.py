#20141106

#A script for producing a SED for freefree emission in W40

from numpy import arange,array,ones, loadtxt#,random,linalg
from pylab import plot,show, xlabel, ylabel, errorbar, log10
from scipy import stats
import matplotlib.pyplot as plt

def line(i,x):
    intercept = data[2][1] - (alpha[i]*data[2][0])
    return alpha[i]*x+intercept

data = loadtxt('data/OS2a_data.txt')

l = len(data[:,0])

x = data[:,0]
erry = data[:,2]
y = data[:,1]

alpha = [0.1,0.6,1.0,1.50] #OS2a

#alpha = [0.21,0.64,0.75,0.82] #OS1a
#alpha = [0.79,1.22,1.32,1.40] #OS2a
#alpha = [0.37,0.80,0.90,0.98] #VLA3a
#alpha = [0.44,0.87,0.97,1.05] #VLA3b
#alpha = [0.40,0.83,0.93,1.01] #VLA3c
#alpha = [0.62,0.86,0.97,1.04] #VLA3d

#print 'intercept', intercept
#print 'slope', slope

LineA = line(0,x) #10%
LineB = line(1,x) #50%
LineC = line(2,x) #75%
LineD = line(3,x) #100%

plot(x,LineA,'r--',label='10 a=%s'%(alpha[0]))
plot(x,LineB,'g--',label='50 a=%s'%(alpha[1]))
plot(x,LineC,'b--',label='75 a=%s'%(alpha[2]))
plot(x,LineD,'k--',label='100 a=%s'%(alpha[3]))
errorbar(x,y,xerr=0,yerr=erry/2,fmt='b+')
plt.legend(loc=2)
xlabel('log Frequency (GHz)')
ylabel('log Flux density (Jy)')
show()
