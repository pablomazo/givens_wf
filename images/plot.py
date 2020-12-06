import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('wf.dat')

data[:,0] *= 1e10 # Convert to A
data[:,1:3] *= 6241506479963200000 # Convert to eV

plt.xlim(0.8,3)
plt.ylim(0,4)
plt.xlabel('$d/\AA$')
plt.ylabel('$E$/eV')
plt.plot(data[:,0], data[:,1], c='black')
for i in range(3,10):
    plt.plot(data[:,0], (data[:,i]/1e6 + data[i,2]), c='blue')
plt.show()
