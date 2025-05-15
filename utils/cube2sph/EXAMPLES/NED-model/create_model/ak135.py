import numpy as np 
import matplotlib.pyplot as plt 
from tqdm import tqdm

def get_gll(NGLL):
    from libgll import gauss_legendre_lobatto,lagrange_poly
    xgll,wgll = gauss_legendre_lobatto(NGLL)
    hprime = np.zeros((NGLL,NGLL))

    for i in range(NGLL):
        _,hprime[i,:] = lagrange_poly(xgll[i],xgll)

    return xgll,wgll,hprime

def compute_jacobian(skel:np.ndarray,xi_loc):
    from libgll import lagrange_poly
    control_loc = np.array([-1,1.])

    poly_termx,poly_derix = lagrange_poly(xi_loc,control_loc)

    # loop every point
    xloc_sum = np.sum(skel * poly_termx)
    dx_dxi = np.sum(skel * poly_derix)
    jaco = dx_dxi 
    xix = 1. / dx_dxi 

    return xloc_sum,jaco,xix

def ak135(r0, param, idom):

    r = r0 / 1000.
    x_ak = r / 6371. #! normalized


    if idom==1 :
        ro_ak = 2.72
        vp_ak = 5.8
        vs_ak = 3.46
        Qmu = 600.0
        Qkappa = 57827.0
    elif idom==2:
        ro_ak = 2.92
        vp_ak = 6.5
        vs_ak = 3.85
        Qmu = 600.0
        Qkappa = 57827.0
    elif idom==3: 
        #! moho -> 210
        ro_ak =  7.1576 - 3.859  * x_ak
        vp_ak = 17.4734 - 9.5332 * x_ak
        vs_ak =  5.8556 - 1.3825 * x_ak
        Qmu = 600.0
        Qkappa = 57827.0
    elif idom==4:
        #! 210 -> 410
        ro_ak =  7.1594 -  3.8608 * x_ak
        vp_ak = 30.7877 - 23.2542 * x_ak
        vs_ak = 15.2181 - 11.0601 * x_ak
        Qmu = 143.0
        Qkappa = 57827.0
    elif idom==5:
        #! 410 -> 660
        ro_ak = 11.1204 -  7.8713 * x_ak
        vp_ak = 29.389  - 21.4066 * x_ak
        vs_ak = 17.7173 - 13.5065 * x_ak
        Qmu = 143.0
        Qkappa = 57827.0
    elif idom==6:
        #! 660 -> D''
        ro_ak =  6.8294 - 1.7227  * x_ak -  1.1064 * x_ak**2 -  0.034409 * x_ak**3
        vp_ak = 26.8598 - 48.9644 * x_ak + 63.7326 * x_ak**2 - 32.4155   * x_ak**3
        vs_ak = 18.0019 - 43.6346 * x_ak + 60.4205 * x_ak**2 - 29.689    * x_ak**3
        Qmu = 312.0
        Qkappa = 57827.0
    elif idom==7:
        #! D'' -> CMB
        ro_ak = -65.8145 + 386.221  * x_ak - 691.6551 *x_ak**2 + 409.6742 * x_ak**3
        vp_ak =   3.4872 + 55.1872  * x_ak -  99.0089 *x_ak**2 +  58.7141 * x_ak**3
        vs_ak = -22.9553 + 164.0287 * x_ak - 294.2766 *x_ak**2 + 174.5113 * x_ak**3
        Qmu = 312.0
        Qkappa = 57827.0
    elif idom==8:
        #! CMB -> ICB
        ro_ak = 12.592  - 1.778  * x_ak - 1.6964 * x_ak**2 -  7.3524 * x_ak**3
        vp_ak = 10.7738 - 2.4831 * x_ak + 3.2584 * x_ak**2 - 14.9171 * x_ak**3
        vs_ak = 0
        Qmu = 0.0
        Qkappa = 57827.0
    elif idom==9:
        #! inner core
        ro_ak = 13.0122 - 0.0011863 *x_ak - 8.4449 * x_ak**2
        vp_ak = 11.2641 - 0.090247  *x_ak - 5.7431 * x_ak**2
        vs_ak =  3.6677 + 0.0049932 *x_ak - 4.4808 * x_ak**2
        Qmu = 84.6
        Qkappa = 1327.7

    if param=='rho':
        out = ro_ak * 1000.
    elif param=='v_p':
        out = vp_ak * 1000.
    elif param=='v_s':
        out = vs_ak * 1000.
    elif param=='vpv':
        out = vp_ak * 1000.
    elif param=='vsv':
        out = vs_ak * 1000.
    elif param=='vph':
        out = vp_ak * 1000.
    elif param=='vsh':
        out = vs_ak * 1000.
    elif param=='eta':
        out = 1.
    elif param=='Qmu':
        out = Qmu
    elif param=='Qka':
        out = Qkappa
    else:
        out =None 

    return out 

def main():
    # set parameters
    sigma = 5

    # mesh
    NGLL = 5 
    ak135_discon = np.array([0,20,35,210,410.,660.,2740.,2889]) # up to 2889km
    nspec_discon = np.array([20,20,60,80,60,100,30])
    nspec = np.sum(nspec_discon)
    skel = np.zeros((nspec + 1))
    n = 0
    for id in range(len(nspec_discon)):
        xmin = ak135_discon[id]
        xmax = ak135_discon[id+1]
        dx = (xmax - xmin) / nspec_discon[id]
        for i in range(nspec_discon[id]):
            skel[n] = xmin + i * dx 
            n = n + 1
    skel[n] = ak135_discon[-1]

    # gll
    xgll,wxgll,hprime = get_gll(NGLL)

    # set jacobian
    jaco = np.zeros((nspec,NGLL))
    xix = jaco * 1. 
    xx = np.zeros((nspec,NGLL))
    for ispec in range(nspec):
        for igll in range(NGLL):
            xx[ispec,igll],jaco[ispec,igll],xix[ispec,igll] = compute_jacobian(skel[ispec:ispec+2],xgll[igll])

    # ibool 
    ibool = np.zeros((nspec,NGLL),dtype=int)
    num = 0
    for i in range(nspec):
        for igll in range(NGLL):
            ibool[i,igll] = num
            if igll < NGLL-1:
                num += 1
    nglob = np.max(ibool) + 1

    # rmassinv
    rmassinv = np.zeros((nglob))
    for ispec in range(nspec):
        mass = wxgll * jaco[ispec,:]
        for igll in range(NGLL):
            iglob = ibool[ispec,igll]
            rmassinv[iglob] += mass[igll]
    rmassinv = 1. / rmassinv

    # prepare C tensor
    dmin = np.unique(np.diff(xx.flatten()))[1]
    dt = 1.0 / 12 * dmin * dmin 
    nt = int(sigma**2 / ( 2 * dt)) + 1
    c = dt

    # parameters
    vp = np.zeros((nspec,NGLL))
    vs = np.zeros((nspec,NGLL))
    rho = np.zeros((nspec,NGLL))
    
    print('smoothing begin ...')
    n = 0
    for id in range(len(nspec_discon)):
        for i in range(nspec_discon[id]):
            for igll in range(NGLL):
                r = 1000 * (6371 - xx[n,igll])
                vp[n,igll] = ak135(r,'v_p',id + 1)
                vs[n,igll] = ak135(r,'v_s',id + 1)
                rho[n,igll] = ak135(r,'rho',id + 1)
            n = n + 1
    
    displ_vp = np.zeros((nglob))
    displ_vs = np.zeros((nglob))
    displ_rho = np.zeros((nglob))
    counter = np.zeros((nglob),dtype=int)
    xstore = displ_vp * 0.
    for ispec in range(nspec):
        for igll in range(NGLL):
            iglob = ibool[ispec,igll]
            displ_vp[iglob] += vp[ispec,igll]
            displ_vs[iglob] += vs[ispec,igll]
            displ_rho[iglob] += rho[ispec,igll]
            counter[iglob] += 1
            xstore[iglob] = xx[ispec,igll]
    displ_vp /= counter
    displ_vs /= counter
    displ_rho /= counter

    # smoothing
    for it in tqdm(range(nt)):
        accel_vp = np.zeros((nglob))
        accel_rho = np.zeros((nglob))
        accel_vs = np.zeros((nglob))
        for ispec in range(nspec):
            # u1 = np.zeros((NGLL))
            # u2 = np.zeros((NGLL))
            # u3 = np.zeros((NGLL))
            idx = ibool[ispec,:]
            u1 = displ_vp[idx]
            u2 = displ_vs[idx]
            u3 = displ_rho[idx]
            # for igll in range(NGLL):
            #     u1[igll] = displ_vp[ibool[ispec,igll]]
            #     u2[igll] = displ_vs[ibool[ispec,igll]]
            #     u3[igll] = displ_rho[ibool[ispec,igll]]
        
            s1 = c * np.dot(hprime,u1) * xix[ispec,:]
            s2 = c * np.dot(hprime,u2) * xix[ispec,:]
            s3 = c * np.dot(hprime,u3) * xix[ispec,:]

            f1 = np.dot(s1*jaco[ispec,:] * xix[ispec,:] * wxgll,hprime)
            f2 = np.dot(s2*jaco[ispec,:] * xix[ispec,:] * wxgll,hprime)
            f3 = np.dot(s3*jaco[ispec,:] * xix[ispec,:] * wxgll,hprime)
            accel_vp[idx] -= f1 
            accel_vs[idx] -= f2
            accel_rho[idx] -= f3

            # for igll in range(NGLL):
            #     iglob = ibool[ispec,igll]
            #     accel_vp[iglob] -= f1[igll]
            #     accel_vs[iglob] -= f2[igll]
            #     accel_rho[iglob] -= f3[igll]
        accel_vp *= rmassinv 
        accel_vs *= rmassinv 
        accel_rho *= rmassinv 
        displ_rho += accel_rho
        displ_vp += accel_vp
        displ_vs += accel_vs

    # smoothing
    f = open("ak135.smooth.bm","w")
    f.write("NAME ak135_smooth\n")
    f.write("ANELASTIC       T\n")
    f.write("ANISOTROPIC     T\n")
    f.write('UNITS           m\n')
    f.write("COLUMNS       radius      rho      vpv      vsv      qka      qmu      vph      vsh      eta\n")

    for iglob in range(nglob):
        vp0 = displ_vp[iglob]
        vs0 = displ_vs[iglob]
        rho0 = displ_rho[iglob]
        f.write("\t %f %f %f %f %f %f %f %f %f\n"%(6371000 - 1000*xstore[iglob],rho0,vp0,vs0,
                                                9999.,9999.,vp0,vs0,1.))
    f.write("#          Discontinuity   1, depth:    2889.00 km\n")
    n = 150
    z = np.linspace(3482,1217,n)
    for i in range(len(z)):
        z0 = z[i]
        r = 1000 * z0
        vp0 = ak135(r,'v_p',8)
        vs0= ak135(r,'v_s',8)
        rho0 = ak135(r,'rho',8)
        f.write("\t %f %f %f %f %f %f %f %f %f\n"%(r,rho0,vp0,vs0,
                                                9999.,9999.,vp0,vs0,1.))
    
    f.write("#          Discontinuity   2, depth:    5154.00 km\n")
    n = 150
    z = np.linspace(1217,0,n)
    for i in range(len(z)):
        z0 = z[i]
        r = 1000 * z0
        vp0 = ak135(r,'v_p',9)
        vs0= ak135(r,'v_s',9)
        rho0 = ak135(r,'rho',9)
        f.write("\t %f %f %f %f %f %f %f %f %f\n"%(r,rho0,vp0,vs0,
                                                9999.,9999.,vp0,vs0,1.))
    # write inner
    f.close()

    # save smoothed model
    f = open("ak135.txt","w")
    for i in range(len(displ_rho)):
        f.write("%f %f %f %f\n" %(-xstore[i],displ_rho[i],displ_vs[i],displ_vp[i]))
    f.close()

    plt.figure(1,figsize=(6,14))
    plt.plot(displ_vs/1000,-xstore)
    plt.plot(vs.flatten()/1000,-xx.flatten())
    plt.savefig("smooth.jpg")
main()