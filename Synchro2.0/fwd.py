# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 19:17:09 2020

@author: zhuravlev.sd
"""


import numpy as np
import math, cmath


def fwMT(freq,model, flag):
    
    for i in range (0, len(model)):
        if model[i]<1:
            apparentResistivity=np.zeros((np.shape(freq)))
            return apparentResistivity
    
    resistivities=model[0:int(len(model)/2)+1]
    thicknesses=model[int(len(model)/2)+1:]
    
    print('MT')
    print('r = ', resistivities)
    print('h = ', thicknesses)
    
    mu = 4*math.pi*1E-7
    n = len(resistivities)
    apparentResistivity=np.zeros((np.shape(freq)))
    phase=np.zeros((np.shape(freq)))
    k=0
    for frequency in freq:   
        w =  2*math.pi*frequency       
        impedances = list(range(n))
        #compute basement impedance
        impedances[n-1] = cmath.sqrt(w*mu*resistivities[n-1]*1j)


   
        for j in range(n-2,-1,-1):
            resistivity = resistivities[j]
            thickness = thicknesses[j]
            
    
            dj = cmath.sqrt((w * mu * (1.0/resistivity))*1j)
            wj = dj * resistivity
            
            ej = cmath.exp(-2*thickness*dj)                     
            belowImpedance = impedances[j + 1]
            rj = (wj - belowImpedance)/(wj + belowImpedance)
            re = rj*ej
            Zj = wj * ((1 - re)/(1 + re))
            impedances[j] = Zj    
        Z = impedances[0]
        absZ = abs(Z);
        apparentResistivity[k] = (absZ * absZ)/(mu * w)
        phase[k] = math.atan2(Z.imag, Z.real)
        k=k+1
    if flag ==0:
        return apparentResistivity
    elif flag ==1:
        return apparentResistivity, phase



#_____________________________generate syntetic data___________________________
f = np.logspace(np.log10(10e-5), np.log10(10e+4), 100)
r = [30,25000,10,5000,2.5,1000]
h = [100,100,100,250,500]
mod=r+h


MTDATA_r =fwMT(f,mod,0)

data=np.zeros((len(MTDATA_r), 2))
data[:,0]=f
data[:,1]=MTDATA_r

np.save('mt_data',data)

    