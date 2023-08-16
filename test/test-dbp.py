
import numpy as np
import matplotlib.pyplot as plt

d   = np.arange(1.8, 2.2, 0.01)

#dbp = np.exp(-15.0*(d-1.8)**2)

dbp = 1.0/(1+np.exp(20.0*(d-2.0)))

plt.figure()
plt.plot(d,dbp)
plt.show()


