#20141106

#A script for producing a SED for freefree emission in W40

from numpy import arange,array,ones, loadtxt#,random,linalg
from pylab import plot,show, xlabel, ylabel, errorbar, log10
from scipy import stats
import matplotlib.pyplot as plt

def line(i,x):
    intercept = data[2][1] - (alpha[i]*data[2][0])
    return alpha[i]*x+intercept

def lineS2(i,x):
    #free-free SED
    #intercept = data[2][1] - (alpha[i]*data[2][0])
    #DUST
    intercept = data[1][1] - (alpha[i]*data[1][0])
    return alpha[i]*x+intercept

def alphaf():
    #calc 100% grad
    #return (data[1][1]-data[2][1])/(data[1][0]-data[2][0])
    #DUST
    return (data[0][1]-data[1][1])/(data[0][0]-data[1][0])

data = loadtxt('data/W40_data.txt')
#data = loadtxt('data/OS1a_data.txt')

l = len(data[:,0])

x = data[:,0]
erry = data[:,2]
y = data[:,1]

alpha = [-0.1,0.6,1.0] 
dust = alphaf()
alpha.append(dust)

print data[0][1],data[1][1]

#print 'intercept', intercept
#print 'slope', slope

LineA = line(0,x) #-01
LineB = line(1,x) #0.6
LineC = line(2,x) #1.0
LineD = lineS2(3,x) #dust

plot(x,LineA,'r-',label=r'$\alpha_{\mathrm{ff}}$ = %s'%(alpha[0]))
#plot(x,LineB,'g-',label=r'$\alpha_{\mathrm{ff}}$ = %s'%(alpha[1]))
#plot(x,LineC,'b--',label=r'$\alpha_{\mathrm{ff}}$ = %s'%(alpha[2]))
#plot(x,LineD,'k--',label=r'$\alpha_{\mathrm{ff}}$ = %s'%(round(alpha[3],1))+'\n(100%)')
plot(x,LineD,'k--',label=r'$\alpha_{\mathrm{dust}}$ = %s'%(round(alpha[3],1)))

errorbar(x,y,xerr=0,yerr=erry,fmt='k+')
plt.legend(loc=2)
xlabel('log Frequency (GHz)')
ylabel('log Flux density (Jy)')
show()
