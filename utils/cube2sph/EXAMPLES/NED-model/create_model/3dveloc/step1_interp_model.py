import numpy as np 
import os 
from scipy.interpolate import griddata,interp1d
import matplotlib.pyplot as plt

# set parameters  begin

# module load 
module_cmd = "module load intel openmpi hdf5 netcdf"

# cube2sph region
xmin=-270000
xmax=270000
ymin=-300000
ymax=300000
zmin=-220000
zmax=0.

# cube2sph region without pml + half elements (4.5 elements)
xmin1=-249750.0
xmax1=249750.0
ymin1=-278906.25
ymax1=278906.25
zmin1=-185000

# smoothing parameter in m
SIGMA_H = 8000
SIGMA_V = 8000

# cube2sph transformation params, from step0_plot_region.sh
cube2sph_param="44.3750000000000       -122.250000000000    0.000000000000000E+00"

# out file
outfile = "tomography_model.3d.xyz"

###### STOP HERE !!!!!! #######################

# functions 

def find_loc(x0,x):
    return np.argmin(np.abs(x0 - x))

def find_taperloc(x0,x1,x):
    ix1 = np.argmin(abs(x0 - x))
    if x0 > x[ix1]:
        ix1 = ix1 + 1
    ix2 = ix1 + 4
    ix4 = np.argmin(abs(x1 - x))
    if x1 < x[ix4]:
        ix4 = ix4 - 1
    ix3 = ix4 - 4

    return ix1,ix2,ix3,ix4

def costaper(n):
    return  0.5 * (1 + np.cos(np.linspace(np.pi, 2 * np.pi, n)))

def add_taper(n,ix1,ix2,ix3,ix4):
    taper = np.zeros((n))
    taper[ix1:ix2+1] = costaper(ix2 - ix1 + 1)
    taper[ix3:ix4 + 1] = costaper(ix4-ix3 + 1)[::-1]
    taper[ix2+1:ix3] = 1.

    return taper 

def get_ak135(z):
    # load velocity model
    data = np.loadtxt("ak135.txt")
    data[:,0] *= 1000

    # create interpolant
    vp_int = interp1d(data[:,0],data[:,1])
    vs_int = interp1d(data[:,0],data[:,2])
    rho_int = interp1d(data[:,0],data[:,3])

    # do interpolate
    nz = len(z)
    vp = np.zeros((nz))
    vs = np.zeros((nz))
    rho = np.zeros((nz))
    for i in range(nz):
        vp[i] = vp_int(z[i])
        vs[i] = vs_int(z[i])
        rho[i] = rho_int(z[i])

    return vp,vs,rho 

def gauss_smoothing(sigx,sigy,sigz,data,nz,ny,nx):
    from scipy.ndimage import gaussian_filter
    newdata = data.reshape(nz,ny,nx)
    print("min/max before smoothing: ",newdata.min(),newdata.max())

    newdata = gaussian_filter(newdata,[sigz,sigy,sigx])
    print("min/max after smoothing: ",newdata.min(),newdata.max())

    return newdata.flatten()

def main():

    # get point cloud in top of cube 
    nx = 151
    ny = 151
    nz = 126
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)
    z = np.linspace(zmin,zmax,nz)

    # write file
    fio = open("model_cube.in","w")
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                fio.write("%f %f %f\n" % (y[iy],x[ix],z[iz]))
    fio.close()

    # cube2latlon
    print("cube2sph transformation ... ")
    cmd = module_cmd + ";" + f"../../bin/cube2latlon model_cube.in model_cube.out {cube2sph_param} "
    os.system(cmd)

    # interpolate
    print("interpolating velocity by using model3d.txt ...")

    # load coordinates
    cords = np.loadtxt("model_cube.out")[:,:3]
    idx = cords[:,1] > 180
    cords[idx,1] = cords[idx,1] - 360
    tmp = cords[:,0] * 1.
    cords[:,0] = cords[:,1] * 1.
    cords[:,1] = tmp * 1.
    cords[:,2] = (6371000 - cords[:,2]) / 1000 # depth

    # load ak135 model
    vp_1d,vs_1d,rho_1d = get_ak135(z)

    # load external 3d model
    model = np.loadtxt("model3d.txt")
    vp_3d = griddata(model[:,:3],model[:,3],cords,rescale=True) * 1000
    vs_3d = griddata(model[:,:3],model[:,4],cords,rescale=True) * 1000
    rho_3d = griddata(model[:,:3],model[:,5],cords,rescale=True) * 1000
    idx = np.isnan(vp_3d)
    vp_3d[idx] = griddata(model[:,:3],model[:,3],cords[idx,:],rescale=True,method='nearest') * 1000
    vs_3d[idx] = griddata(model[:,:3],model[:,4],cords[idx,:],rescale=True,method='nearest') * 1000
    rho_3d[idx] = griddata(model[:,:3],model[:,5],cords[idx,:],rescale=True,method='nearest') * 1000

    # smoothing
    print("smoothing model ...")
    sx = SIGMA_H / (xmax - xmin) / 2 * nx
    sy = SIGMA_H / (ymax - ymin) / 2 * ny 
    sz = SIGMA_V / (zmax - zmin) / 2 * nz
    vp_3d = gauss_smoothing(sx,sy,sz,vp_3d,nz,ny,nx)
    vs_3d = gauss_smoothing(sx,sy,sz,vs_3d,nz,ny,nx)
    rho_3d = gauss_smoothing(sx,sy,sz,rho_3d,nz,ny,nx)

    print("taper to ak135 at PML ...")

    # get taper loc
    ix1,ix2,ix3,ix4 = find_taperloc(xmin1,xmax1,x)
    iy1,iy2,iy3,iy4 = find_taperloc(ymin1,ymax1,y)
    iz1,iz2,_,_ = find_taperloc(zmin1,zmax,z)
    taperx = add_taper(nx,ix1,ix2,ix3,ix4)
    tapery = add_taper(ny,iy1,iy2,iy3,iy4)
    taperz = np.zeros((nz))
    taperz[iz1:iz2+1] = costaper(iz2 - iz1 + 1)
    taperz[iz2+1:] = 1.

    #plot
    plt.figure(1,figsize=(14,6))
    plt.subplot(131)
    plt.plot(x,taperx)
    plt.axvline(xmin1)
    plt.axvline(xmax1)
    plt.subplot(132)
    plt.plot(y,tapery)
    plt.axvline(ymin1)
    plt.axvline(ymax1)   
    plt.subplot(133)
    plt.plot(z,taperz)
    plt.axvline(zmin1)   
    plt.savefig("taper.jpg")

    # taper to ak135 at PML
    for iz in range(nz):
        vs0 = vs_1d[iz]
        vp0 = vp_1d[iz]
        rh0 = rho_1d[iz]
        for iy in range(ny):
            for ix in range(nx):
                tp = taperx[ix] * tapery[iy] * taperz[iz]
                idx = iz * ny * nx + iy * nx + ix
                vs = vs0 + (vs_3d[idx] - vs0) * tp 
                rho = rh0 + (rho_3d[idx] - rh0) * tp 
                vp = vp0 + (vp_3d[idx] - vp0) * tp 

                # copy to 3d model
                vs_3d[idx],vp_3d[idx],rho_3d[idx] = vs,vp,rho

    f = open(outfile,"w")
    f.write("%f %f %f %f %f %f\n" %(xmin,ymin,zmin,xmax,ymax,zmax))
    f.write("%f %f %f\n" %(x[1]-x[0],y[1]-y[0],z[1]-z[0]))
    f.write("%d %d %d\n"%(nx,ny,nz))
    f.write("%f %f %f %f %f %f\n" %(vp_3d.min()*0.99,vp_3d.max()*1.01,  \
                                    vs_3d.min()*0.99,vs_3d.max()*1.01,
                                    rho_3d.min()*0.99,rho_3d.max() * 1.01))
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                idx = iz * ny * nx + iy * nx + ix
                f.write("%f %f %f " %(x[ix],y[iy],z[iz]))
                f.write("%f %f %f\n" %(vp_3d[idx],vs_3d[idx],rho_3d[idx]))
    f.close()

    # clean files
    os.remove("model_cube.in")
    os.remove("model_cube.out")


if __name__ == "__main__":
    main()

