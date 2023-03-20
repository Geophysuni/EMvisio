# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 17:13:35 2023

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
from forwards import fwCS, fwMT
# from inverses_new import MTinv, CSinv, Jointinv #, Jointinv2
import numpy.random as rnd
from scipy.spatial import distance
from scipy.optimize import minimize as smin
from scipy.signal import savgol_filter
import empymod
import math
import cmath

f = np.logspace(np.log10(5*10e-5), np.log10(10e+3), 120)
r = np.array([30,25000,10,5000,2.5,1000])
h = np.array([100,100,100,250,500])

# mod=r+h
# CSDATA=fwCS(f,mod)

class Model:
    def __init__(self, thickness, resistivity):
        self.h = thickness
        self.r = resistivity
        
class FProblem:
    def __init__(self, probType, modObj):
        self.type = probType
        self.model = modObj
        self.params = np.nan
        self.measureRange = np.nan
        self.result = np.nan
    
    def setMeasureRange(self, mr):
        self.measureRange = mr
    
    def setParams(self, params):
        self.params = params
    
    def calculate(self):
        if self.type == 'cs':
            # print (self.type)

            resistivities=self.model.r
            thicknesses=self.model.h
            
            depth=np.zeros((np.shape(resistivities)))
            depth[0]=0
            for i in range(1,len(thicknesses)+1):
                depth[i]=depth[i-1]+thicknesses[i-1]
                
            res=np.zeros((len(resistivities)+1))
            
            res[0]=2e14
            res[1:]=resistivities
            
            freq = self.measureRange
            
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
            
            
            r=np.sqrt((model['src'][0]-model['rec'][0])**2+(model['src'][1]-model['rec'][1])**2
                      +(model['src'][2]-model['rec'][2])**2)
            
            CS_rho=res_0*2*math.pi*(r**3)/(3-2)
            
            self.result = abs(CS_rho)
            
        if self.type == 'mt':
            # print (self.type)
            
            resistivities=self.model.r
            thicknesses=self.model.h
        
            freq = self.measureRange
            
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
                
            self.result = apparentResistivity
        
# model = Model(h,r)
# prb = FProblem('mt', model)
# prb.setMeasureRange(f)

plt.figure()

for i in range(10):
    r = r*(np.random.rand()+0.5)
    h = h*(np.random.rand()+0.5)
    model = Model(h,r)
    prb = FProblem('cs', model)
    prb.setMeasureRange(f)
    prb.calculate()
    plt.loglog(f,prb.result)
    
# m_model = 
