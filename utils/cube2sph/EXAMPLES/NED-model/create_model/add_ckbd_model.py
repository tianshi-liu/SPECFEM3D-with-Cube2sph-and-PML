import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm
import matplotlib.pyplot as plt 

def find_loc(x,x0):
    n = len(x)
    dx = x[1] - x[0]

    iloc = int((x0 - x[0]) / dx)

    return iloc

# study region
xmin,xmax = -180000, 180000
ymin,ymax = -180000, 180000
zmin,zmax = -60000,0.

nx,ny,nz = 121,121,51
x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)
z = np.linspace(zmin,zmax,nz)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
print("dx dy dz = ",dx,dy,dz)

# load 1-D model
data = np.loadtxt("ak135.txt")
data[:,0] *= 1000

# create interpolant
vp_int = interp1d(data[:,0],data[:,1])
vs_int = interp1d(data[:,0],data[:,2])
rho_int = interp1d(data[:,0],data[:,3])

# model
vs = np.zeros((nz))
rho = vs * 1.
vp = vs + 0.
for iz in range(nz):
    z0 = z[iz]
    vp0 = vp_int(z0)
    vs0 = vs_int(z0)
    rho0 = rho_int(z0)
    vs[iz] = vs0 
    vp[iz] = vp0 
    rho[iz] = rho0 

# start end index for 
ix0 = find_loc(x,-75000) + 2
ix1 = find_loc(x,75000) - 2
iy0 = find_loc(y,-75000) + 2
iy1 = find_loc(y,75000) - 2
iz0 = find_loc(z,-30000) + 2
iz1 = nz - 2
print(ix0,ix1,iy0,iy1,iz0,iz1)

# size of checkerboard
cksize_x = 35000
cksize_z = 15000
nx0 = int(cksize_x / dx * 2)
ny0 = int(cksize_x / dy * 2)
nz0 = int(cksize_z / dz)

out = np.zeros((ny,nx))
out1 = np.zeros((nz,nx))

# print(x[ix0],y[iy0],z[iz0])
# print(x[ix1],y[iy1],z[iz1])

f = open("tomography_model.xyz","w")
f.write("%f %f %f %f %f %f\n" %(xmin,ymin,zmin,xmax,ymax,zmax))
f.write("%f %f %f\n" %(dx,dy,dz))
f.write("%d %d %d\n"%(nx,ny,nz))
f.write("%f %f %f %f %f %f\n" %(vp.min()*0.9,vp.max()*1.1,vs.min()*0.9,vs.max()*1.1,rho.min()*0.9,rho.max()*1.1))
for iz in tqdm(range(nz)):
    vp0 = vp[iz]
    vs0 = vs[iz]
    rho0 = rho[iz]
    for iy in range(ny):
        for ix in range(nx):
            vp1 = vp0 * 1. 
            vs1 = vs0 * 1.
            rho1 = rho0 * 1.
            if ix >= ix0 and ix <= ix1 and iy >= iy0 and iy <=iy1 and iz >= iz0 and iz<=iz1:
                dv = 0.1 * np.sin(np.pi * (ix- ix0)/ nx0) *   \
                    np.sin(np.pi * (iy - iy0) / ny0) *   \
                    np.sin(np.pi * (iz - iz0) / nz0)
                vp1 = vp0 * (1 + dv)
                vs1 = vs0 * (1 + dv)
                rho1 = rho0 * (1 + dv)
                out[iy,ix] = dv

                out1[iz,ix] = dv

            f.write("%f %f %f %f %f %f\n" %(x[ix],y[iy],z[iz],vp1,vs1,rho1))
f.close()

plt.figure(1,figsize=(12,3))
plt.subplot(121)
plt.contourf(x,y,out)
plt.subplot(122)
plt.contourf(x,z,out1)
plt.savefig("checkboard.jpg")

