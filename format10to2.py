# -*- coding: utf-8 -*-
import math
import numpy as np
f=open('input.txt','w')
g=open("inputtest.txt").readlines()
y=np.zeros(len(g))
z=np.zeros(len(g))
print(y)
#y=np.array([0.66379454, -0.59497408, -2.15299592,  2.84935889,  6.53625405,  2.87575803])
#上位八桁整数、下位24桁少数
#print(g[0])
for j in range (len(g)):
    y[j]=float(g[j])
#print (y)
for i in range(len(g)):
    if y[i]>=0.0:
        #print(format(int(y[i]*(2**24)),'032b'))
        f.write(format(int(y[i]*(2**24)),'032b').rstrip('\r\n')+' ')
    else :
        #print(format(int(2**32+y[i]*(2**24)),'032b'))
        f.write(format(int(2**32+y[i]*(2**24)),'032b').rstrip('\r\n')+' ')


    


