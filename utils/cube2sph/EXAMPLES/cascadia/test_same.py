import numpy as np 
from scipy.io import FortranFile

f1 = FortranFile("KERNEL/proc000000_vp_smooth.bin","r")
a_cpu = f1.read_reals('f4')
f1.close()

f1 = FortranFile("KERNEL_GPU/proc000000_vp_smooth.bin","r")
a_gpu = f1.read_reals('f4')
f1.close()

print(np.allclose(a_gpu,a_cpu))
print(np.max(abs(a_gpu-a_cpu)))

idx = np.argsort(abs(a_gpu - a_cpu))[::-1]
print(a_gpu[idx[:10]])
print(a_cpu[idx[:10]])
