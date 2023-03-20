# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 10:49:01 2020

@author: zhuravlev.sd
"""

import empymod
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

import time

##_____________________________________________________________________________
def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)

##_____________________________________________________________________________

def fwCS(freq,mdl):
    for i in range (0, len(mdl)):
        if mdl[i]<1:
            CS_rho=np.zeros((np.shape(freq)))
            return abs(CS_rho)
    
    resistivities=mdl[0:int(len(mdl)/2)+1]
    thicknesses=mdl[int(len(mdl)/2)+1:]

#    print('CS')
#    print('r = ', resistivities)
#    print('h = ', thicknesses)    
    
    
    depth=np.zeros((np.shape(resistivities)))
    depth[0]=0
    for i in range(1,len(thicknesses)+1):
        depth[i]=depth[i-1]+thicknesses[i-1]
        
    res=np.zeros((len(resistivities)+1))
    
    res[0]=2e14
    res[1:]=resistivities
    
    # Model
    model = {
        'src': [0, 0, 0.001],           # Src at origin, slightly below interface
        'rec': [6000, 0, 0.0001],       # 6 km off., in-line, slightly bel. interf.
        'depth': depth,       # Target of 100 m at 2 km depth
        'res': res,     # Res: [air, overb., target, half-space]
        'epermH': np.ones((len(res))),         # epermH/epermV: correspond to the default
        'epermV': np.ones((len(res))),         # #              values if not provided
        'freqtime': freq,                  # Times
        'signal': None,                    # 0: Impulse response
        'ab': 11,
        'ftarg': 'key_81_CosSin_2009',  # Choose a shorter filter then the default
        'verb': 1,                      # Verbosity; set to 3 to see all parameters
    }
    
    model['epermH'][0] = 0
    model['epermV'][0] = 0
    
    res_0 = empymod.dipole(**model)    
    CS_rho=np.zeros((np.shape(res_0))) 
    
#    r=np.sqrt((model['src'][0][0]-model['rec'][0][0])**2+(model['src'][1][0]-model['rec'][1][0])**2
#              +(model['src'][2][0]-model['rec'][2][0])**2)
    
    r=np.sqrt((model['src'][0]-model['rec'][0])**2+(model['src'][1]-model['rec'][1])**2
              +(model['src'][2]-model['rec'][2])**2)
    
    CS_rho=res_0*2*math.pi*(r**3)/(3-2)
    return abs(CS_rho)

def fwMT(freq,model, flag):
    
    for i in range (0, len(model)):
        if model[i]<1:
            apparentResistivity=np.zeros((np.shape(freq)))
            return apparentResistivity
    
    resistivities=model[0:int(len(model)/2)+1]
    thicknesses=model[int(len(model)/2)+1:]
    
#    print('MT')
#    print('r = ', resistivities)
#    print('h = ', thicknesses)
    
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

def jointFW(freq,model):
    d1=fwMT(freq,model,0)
    d2=fwCS(freq,model)
    d=np.append(d1,d2)
    return d

