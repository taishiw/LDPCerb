# -*- coding: utf-8 -*-
import math
import numpy as np
f=open('input16.txt','w')
g=open("inputtest.txt").readlines()
h=open('outputalpha16.txt','w')
k=open('outputsoftalpha.txt').readlines()
l=open('outputbeta16.txt','w')
m=open('outputsoftbeta.txt').readlines()
x=np.zeros(len(g))
y=np.zeros(len(k))
z=np.zeros(len(m))
#y=np.array([0.66379454, -0.59497408, -2.15299592,  2.84935889,  6.53625405,  2.87575803])
#上位2桁整数、下位6桁少数
#print(g[0])
for j in range (len(g)):
    x[j]=float(g[j])
#print (y)
for i in range(len(g)):
    if x[i]>=0.0:
        f.write(hex(int(x[i]*(16**6)))+'\n')
    else :
        f.write(hex(int(16**8+x[i]*(16**6)))+'\n')

for j in range (len(k)):
    y[j]=float(k[j])
#print (y)
for i in range(len(k)):
    if y[i]>=0.0:
        h.write(hex(int(y[i]*(16**6)))+'\n')
    else :
        h.write(hex(int(16**8+y[i]*(16**6)))+'\n')

for j in range (len(m)):
    z[j]=float(m[j])
#print (y)
for i in range(len(m)):
    if z[i]>=0.0:
        l.write(hex(int(z[i]*(16**6)))+'\n')
    else :
        l.write(hex(int(16**8+z[i]*(16**6)))+'\n')


    


