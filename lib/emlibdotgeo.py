# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 18:21:21 2023

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio as iio
import math
import cmath
import empymod

# asis_image = iio.imread("oblaka-1.jpg")
# ini_image = np.exp(iio.imread("oblaka-1.jpg")[:, :, 0]/100)*100
# f = np.logspace(np.log10(5*10e-5), np.log10(10e+3), 120)

class Model:
    def __init__(self, mod):
        self.model = mod
        self.shape = np.shape(self.model)
        self.xcrd = np.nan
        self.ycrd = np.nan
        self.dx = np.nan
        self.dy = np.nan
        
    def setSize(self, xrange, yrange):
        self.dx = (xrange[-1]-xrange[0])/self.shape[1]
        self.dy = (yrange[-1]-yrange[0])/self.shape[0]
        
        self.xcrd = np.arange(xrange[0], xrange[-1], self.dx)
        self.ycrd = np.arange(yrange[0], yrange[-1], self.dy)
        
class FProblem:
    def __init__(self, probType):
        self.type = probType
        self.model = np.nan
        self.src = np.nan
        self.rec = np.nan
        self.measureRange = np.nan
        self.ystep = np.nan
        self.measuredValues = np.nan
        
    def setModel(self, modObj):
        self.model = modObj
        self.ystep = self.model.dy
    def setYStep(self, yst_val):
        self.ystep = yst_val
    def setSrc(self, srcCrdList):
        self.src = srcCrdList
    def setRec(self, recCrdList):
        self.rec = recCrdList
    def setMeasurRenge(self, measureRange):
        self.measureRange = measureRange
        
    def get1DModel(self, xcrd):
        xshift = min(self.model.xcrd)
        tmpModX = self.model.xcrd - xshift
        tmpXCrd = xcrd - xshift
        
        # print(xshift)
        # print(tmpModX)
        # print(tmpXCrd)
        # print(xcrd)
        # print(xshift)
        
        xInd = int((tmpXCrd-tmpModX[0])/self.model.dx)
        iniModel = self.model.model[:,xInd]
        
        thickness = np.ones(len(iniModel))*self.ystep
        resistivity = np.zeros(len(thickness)+1)
        
                
        # for i in range(len(thickness)):
            # resistivity[i] = np.mean(iniModel[int(i*self.ystep):int((i+1)*self.ystep)])
            
        resistivity[:-1] = iniModel
        resistivity[-1] = iniModel[-1]
 
        
        # print(xInd)
        
        # print(resistivity)
        
        return thickness, resistivity
    
    def calculate(self):
        
        self.measuredValues = np.zeros((len(self.rec), len(self.measureRange)))
        # mlist = []
        
        for irec, crec in enumerate(self.rec):
            # print(self.rec)
            ch, cr = self.get1DModel(crec)
            # print(crec)
            # print(cr)
            if self.type == 'mt':
                freq = self.measureRange
                mu = 4*math.pi*1E-7
                n = len(cr)
                apparentResistivity=np.zeros((np.shape(freq)))
                phase=np.zeros((np.shape(freq)))
                k=0
                for frequency in freq:   
                    w =  2*math.pi*frequency       
                    impedances = list(range(n))
                    impedances[n-1] = cmath.sqrt(w*mu*cr[n-1]*1j)
                    for j in range(n-2,-1,-1):
                        resistivity = cr[j]
                        thickness = ch[j]
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
                # print(apparentResistivity)
                   
                self.measuredValues[irec,:] = apparentResistivity
                # mlist.append(apparentResistivity)
                
            elif self.type == 'cs':
                
                depth=np.zeros((np.shape(cr)))
                depth[0]=0
                for i in range(1,len(ch)+1):
                    depth[i]=depth[i-1]+ch[i-1]  
                res=np.zeros((len(cr)+1))
                res[0]=2e14
                res[1:]=cr
                            
                # Model
                model = {
                    'src': [0, 0, 0.001],           # Src at origin, slightly below interface
                    # 'src': [-3000,3000, 0,0, 0.001, 0.001], 
                    'rec': [6000, 0, 0.0001],
                    # 'rec': [-1000, 1000, 0,0, 0.001, 0.001], # 6 km off., in-line, slightly bel. interf.
                    'depth': depth,       # Target of 100 m at 2 km depth
                    'res': res,     # Res: [air, overb., target, half-space]
                    'epermH': np.ones((len(res))),         # epermH/epermV: correspond to the default
                    'epermV': np.ones((len(res))),         # #              values if not provided
                    'freqtime': self.measureRange,                  # Times
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
                
                self.measuredValues[irec,:] = CS_rho

        

# model = Model(ini_image)
# model.setSize([0,1200], [0,338])

# prb = FProblem('cs')
# prb.setModel(model)
# prb.setMeasurRenge(f)
# prb.setRec([10])

# # prb.setYStep(1)

# prb.calculate()

# plt.figure()
# for i in prb.measuredValues:
#     plt.loglog(prb.measureRange, i)
