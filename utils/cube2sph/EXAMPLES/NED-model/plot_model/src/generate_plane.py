import numpy as np
import sys

def delsph(lat1,lon1,lat2,lon2):
    R = 6371.0 
    from numpy import sin,cos,sqrt 
    from math import atan2

    colatrad1 = np.deg2rad(90 - lat1)
    colatrad2 = np.deg2rad(90 - lat2)
    lonrad1 = np.deg2rad(lon1)
    lonrad2 = np.deg2rad(lon2) 
    dlat = colatrad2 - colatrad1
    dlon = lonrad2 - lonrad1
    lat1rad = np.pi * 0.5 - colatrad1
    lat2rad =  np.pi * 0.5 - colatrad2
    a = sin(dlat * 0.5 ) * sin(dlat * 0.5 ) + sin(dlon * 0.5)  \
        * sin(dlon * 0.5) * cos(lat1rad) * cos(lat2rad)
    c = 2.0 * atan2(sqrt(a),sqrt(1-a))
    

    return R * c    

def main():
    if len(sys.argv)!=7:
        print("Usage ./this lonmin lonmax latmin latmax n outfile")
        exit(1)
    
    lonmin,lonmax,latmin,latmax = map(lambda x:float(x),sys.argv[1:5])
    n = int(sys.argv[5])
    outfile = sys.argv[6]
    data = np.zeros((n*n,2))

    for i in range(n):
        for j in range(n):
            data[i*n+j,0] = lonmin + (lonmax - lonmin) / (n-1) * i 
            data[i*n+j,1] = latmin + (latmax - latmin) / (n-1) * j 
    np.savetxt(outfile,data,fmt='%f')

main()
