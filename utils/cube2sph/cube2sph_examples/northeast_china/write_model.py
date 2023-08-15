import numpy as np

x0 = -1111949.266
y0 = -555974.633
dx = 2223898.532
dy = 1111949.266
nx = 1
ny = 1

rho_min = 10000.0
vp_min = 10000.0
vs_min = 10000.0
rho_max = 0.0
vp_max = 0.0
vs_max = 0.0

z0 = -1060000.0
dz = 5000.0
nz = 80

A = np.loadtxt("model_ak135")
A[:,0] = -A[:,0]
l1 = 0
l2 = 9

with open("model_bot_660.tomo",'w') as f:
  for k in range(nz+1):
    z = z0+k*dz
    l = l1
    while (l < l2-1):
      if (A[l,0]<=z<=A[l+1,0]): break
      l = l + 1
    rz = (z-A[l,0]) / (A[l+1,0]-A[l,0])
    rho = A[l,1] + rz * (A[l+1,1]-A[l,1])
    rho = rho * 1000.0
    vp = A[l,2] + rz * (A[l+1,2]-A[l,2])
    vp = vp * 1000.0
    vs = A[l,3] + rz * (A[l+1,3]-A[l,3])
    vs = vs * 1000.0
      
    for j in range(ny+1):
      y = y0+j*dy
      for i in range(nx+1):
        x = x0+i*dx
        if (rho > rho_max): rho_max = rho
        if (rho < rho_min): rho_min = rho
        if (vp > vp_max): vp_max = vp
        if (vp < vp_min): vp_min = vp
        if (vs > vs_max): vs_max = vs
        if (vs < vs_min): vs_min = vs
        f.write(f"{x} {y} {z} {vp:.1f} {vs:.1f} {rho:.1f}\n")

with open("model_bot_660.head", "w") as f:
  f.write(f"{x0} {y0} {z0} {x0+nx*dx} {y0+ny*dy} {z0+nz*dz}\n")
  f.write(f"{dx} {dy} {dz}\n")
  f.write(f"{nx+1} {ny+1} {nz+1}\n")
  f.write(f"{vp_min} {vp_max} {vs_min} {vs_max} {rho_min} {rho_max}\n")
