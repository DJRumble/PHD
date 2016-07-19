#Test for producing normal distributions

mean = 2
sigma = 0.1

import random
import numpy as np
import matplotlib.pyplot as plt

i = 0
r = []

while i  in range(5000):
    a = random.normalvariate(mean,sigma)
    r.append(a)
     
    i = i + 1

#print r

H = np.histogram(r,11)

#print H[1]
print H[0]

plt.plot(H[1], H[0])
plt.show()
