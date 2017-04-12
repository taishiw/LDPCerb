# -*- coding: utf-8 -*-
import math
import numpy as np
f=open('outputalpha10.txt','w')
g=open('outputbeta10.txt','w')
h=open("outputalpha.txt").readlines()
l=open("outputbeta.txt").readlines()
y=np.zeros(len(h))
z=np.zeros(len(l))
k=np.zeros(len(h))
m=np.zeros(len(l))
#y=np.array([0.66379454, -0.59497408, -2.15299592,  2.84935889,  6.53625405,  2.87575803])
#上位八桁整数、下位24桁少数
print(l[0])
for j in range (len(h)):
    y[j]=int(h[j],2)
for i in range(len(h)):
    if y[i]<=2**31:
        #print(y[i]*(2**(-24)))
        f.write(str(y[i]*(2**(-24))).rstrip('\r\n')+'\n')
    else :
        #z[i]=format(int(2**32-y[i]),'032b')
        k[i]=int(format(int(2**32-y[i]),'032b'),2)
        #print((-1)*((k[i]*(2**(-24)))))
        f.write(str((-1)*((k[i]*(2**(-24))))).rstrip('\r\n')+'\n')
for j in range (len(l)):
    z[j]=int(l[j],2)
for i in range(len(l)):
    if z[i]<=2**31:
        #print(y[i]*(2**(-24)))
        g.write(str(z[i]*(2**(-24))).rstrip('\r\n')+'\n')
    else :
        #z[i]=format(int(2**32-y[i]),'032b')
        m[i]=int(format(int(2**32-z[i]),'032b'),2)
        #print((-1)*((k[i]*(2**(-24)))))
        g.write(str((-1)*((m[i]*(2**(-24))))).rstrip('\r\n')+'\n')


    


