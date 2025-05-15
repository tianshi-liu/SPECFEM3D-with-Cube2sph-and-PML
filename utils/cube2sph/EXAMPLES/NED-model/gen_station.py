import numpy as np 

x = np.linspace(-1,1,7)
y = x * 1. 

fio = open("station.lst","w")
for i in range(len(x)):
    for j in range(len(x)):
        fio.write("%f %f\n" %(x[i],y[j]))
fio.close()