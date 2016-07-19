#20141106

#A script for producing a SED for freefree emission in W40

from numpy import arange,array,ones, loadtxt#,random,linalg
from pylab import plot,show, xlabel, ylabel, errorbar, log10
from scipy import stats

def line(x):
    return slope*x+intercept

data = loadtxt('free-free_data.txt')

l = len(data[:,0])

x = data[:,1]
erry = data[:,2]*0.1
wav = data[:,0]
y = data[:,2]

slope, intercept, r_value, p_value, slope_std_err = stats.linregress(x,y)

print 'intercept', intercept
print 'slope', slope
print 'standard deviation', slope_std_err
Line = line(x)

v450 = log10(666)
v850 = log10(352)

s450 = slope*v450+intercept
s850 = slope*v850+intercept

vS2 = [v450,v850]
fS2 = [s450,s850]
sS2 = [log10(637),log10(66)]

Y = [Line[0],line(v450),line(v850)]
X = [data[0][1],v450,v850]

plot(X,Y,'r-',x,y,'b+')
plot(vS2,fS2,'b+')
plot(vS2,sS2,'b+')
errorbar(x,y,xerr=0,yerr=erry/2,fmt='b+')

xlabel('log Frequency (GHz)')
ylabel('log Flux density (Jy)')
show()
