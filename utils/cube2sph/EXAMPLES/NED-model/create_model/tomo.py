import numpy as np
from scipy.interpolate import interp1d


# study region
xmin,xmax = -180000, 180000
ymin,ymax = -180000, 180000
zmin,zmax = -60000,0.

nx,ny,nz = 2,2,101
x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)
z = np.linspace(zmin,zmax,nz)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
print("dx dy dz = ",dx,dy,dz)

# load velocity model
data = np.loadtxt("ak135.txt")
data[:,0] *= 1000

# create interpolant
vp_int = interp1d(data[:,0],data[:,1])
vs_int = interp1d(data[:,0],data[:,2])
rho_int = interp1d(data[:,0],data[:,3])

# model
vs = np.zeros((nz,ny,nx))
rho = vs * 1.
vp = vs + 0.
for iz in range(nz):
    z0 = z[iz]
    vp0 = vp_int(z0)
    vs0 = vs_int(z0)
    rho0 = rho_int(z0)
    vs[iz,:,:] = vs0 
    vp[iz,:,:] = vp0 
    rho[iz,:,:] = rho0 

f = open("tomography_model.xyz","w")
f.write("%f %f %f %f %f %f\n" %(xmin,ymin,zmin,xmax,ymax,zmax))
f.write("%f %f %f\n" %(dx,dy,dz))
f.write("%d %d %d\n"%(nx,ny,nz))
f.write("%f %f %f %f %f %f\n" %(vp.min(),vp.max(),vs.min(),vs.max(),rho.min(),rho.max()))
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            f.write("%f %f %f %f %f %f\n" %(x[ix],y[iy],z[iz],vp[iz,iy,ix],vs[iz,iy,ix],rho[iz,iy,ix]))
f.close()

