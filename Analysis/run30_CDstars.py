import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy import stats

def line(intercept,slope,x):
    X = np.array(x)
    y = (slope*X)+intercept
    return y

plt.figure()

bins = [5.62341325,12.11527659,26.10157216,56.23413252,121.15276586,261.0157215, 562.34132519, 1211.52765863]

PRE = [   8.,   75.,  125.,   46.,   11.,    3.,    1.,    0.]
PRO = [  1.,  27.,  68.,  64.,  33.,  19.,   7.,   2.]
TOT = [   9.,  102.,  193.,  110.,   44.,   22.,    8.,    2.]

percent_PRE = []
percent_PRO = []

for i in range(len(TOT)):
    percent_PRE.append(PRE[i]/TOT[i])
    percent_PRO.append(PRO[i]/TOT[i])

print bins
print percent_PRE

plt.semilogx(bins, percent_PRE,'blue',linestyle='none',marker='o',linewidth=1)
plt.semilogx(bins, percent_PRO,'K',linestyle='none',marker='o',linewidth=1)

logbins = np.log(bins)

slope, intercept, r_value1, p_value1, std_err1 = stats.linregress(bins,percent_PRE)
X = [5,3000]
intercept = 1.2
plot(X,line(intercept,slope,X),'b')

slope, intercept, r_value1, p_value1, std_err1 = stats.linregress(bins,percent_PRO)
X = [5,3000]
intercept = -0.2
plot(X,line(intercept,slope,X),'k')

xlim([5,3000])
ylim([0,1])

plt.ylabel('Fraction of clumps in bin')
plt.xlabel('Peak column density (10$^{21}$ H$_{2}$ cm$^{-2}$)')

plt.show()
