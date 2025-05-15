from scipy.io import FortranFile
import numpy as np 

iproc = 100
filename = 'proc%06d_vp.bin' %iproc
f = FortranFile(f'DATABASES_MPI/{filename}',"r")
f1 = FortranFile(f'DATABASES_MPI_REF/{filename}',"r")
v = f.read_reals('f4')
v1 = f1.read_reals('f4')
f.close()
f1.close()

print(np.linalg.norm(v1-v),np.min(v1),np.min(v))
