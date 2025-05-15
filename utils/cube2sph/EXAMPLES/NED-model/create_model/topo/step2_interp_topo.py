import numpy as np 
import os 
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

# set parameters  begin

# module load 
module_cmd = "module load intel openmpi hdf5 netcdf"

# cube2sph region
xmin=-270000
xmax=270000
ymin=-300000
ymax=300000

# cube2sph region without pml
xmin1=-249750.0
xmax1=249750.0
ymin1=-278906.25
ymax1=278906.25

# cube2sph transformation params, from step0_plot_region.sh
cube2sph_param="44.3750000000000       -122.250000000000    0.000000000000000E+00"

# out file
outfile = "interface_top.dat"

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

def main():

    # get point cloud in top of cube 
    nx=151
    ny=151
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)

    # write file
    fio = open("topo.in","w")
    for iy in range(ny):
        for ix in range(nx):
            fio.write("%f %f %f\n" % (y[iy],x[ix],0.))
    fio.close()

    # cube2latlon
    print("cube2sph transformation ... ")
    cmd = module_cmd + ";" + f"../../bin/cube2latlon topo.in topo.out {cube2sph_param} "
    os.system(cmd)

    # interpolate
    print("interpolating topography by using external_topo.txt ...")

    # load coordinates
    cords = np.loadtxt("topo.out")[:,:2]
    idx = cords[:,1] > 180
    cords[idx,1] = cords[idx,1] - 360
    tmp = cords[:,0] * 1.
    cords[:,0] = cords[:,1] * 1.
    cords[:,1] = tmp * 1.

    # load external_topo
    topo = np.loadtxt("external_topo.txt")
    topo_intp = griddata(topo[:,:2],topo[:,2],cords)

    print("taper topography to 0 at PML")

    # get taper loc
    ix1,ix2,ix3,ix4 = find_taperloc(xmin1,xmax1,x)
    iy1,iy2,iy3,iy4 = find_taperloc(ymin1,ymax1,y)
    taperx = add_taper(nx,ix1,ix2,ix3,ix4)
    tapery = add_taper(ny,iy1,iy2,iy3,iy4)

    # plot
    # plt.figure(1,figsize=(14,6))
    # plt.subplot(121)
    # plt.plot(x,taperx)
    # plt.axvline(xmin1)
    # plt.axvline(xmax1)
    # plt.subplot(122)
    # plt.plot(y,tapery)
    # plt.axvline(ymin1)
    # plt.axvline(ymax1)   
    # plt.savefig("taper.jpg")

    # taper topography
    fio = open(outfile,"w")
    for iy in range(ny):
        for ix in range(nx):
            idx = iy * nx + ix
            z = topo_intp[idx] * taperx[ix] * tapery[iy]
            fio.write("%g\n" %(z))
    fio.close()

    print(f"\n\n.true. %d %d %f %f %f %f\n {outfile}" %(nx,ny,x[0],y[0],x[1]-x[0],y[1]-y[0]))

    # write topography in lon/lat
    fio = open("tapered_topo.txt","w")
    for iy in range(ny):
        for ix in range(nx):
            idx = iy * nx + ix
            z = topo_intp[idx] * taperx[ix] * tapery[iy]
            fio.write("%g %g %g\n" %(cords[idx,0],cords[idx,1],z))
    fio.close()

    # clean files
    os.remove("topo.in")
    os.remove("topo.out")


if __name__ == "__main__":
    main()

