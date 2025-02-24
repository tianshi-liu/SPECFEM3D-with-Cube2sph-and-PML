import numpy as np 
import matplotlib.pyplot as plt 
import obspy

name = "TA.K04A.MXZ.semd"
name = "XC.OR012.MXZ.semd"

d = np.loadtxt(f"OUTPUT_FILES/{name}")
d_cpu = np.loadtxt(f"OUTPUT_FILES.cpu/{name}")

# tr = obspy.Trace(d[:,1])
# tr.stats.delta = d[1,0] - d[0,0]
# tr.detrend('linear')
# d[:,1] = tr.data

plt.plot(d[:,0],d[:,1],label='gpu')
plt.plot(d_cpu[:,0],d_cpu[:,1],label='cpu',ls='--')
plt.legend()

plt.savefig("test.jpg")

