# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 20:30:25 2020

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from forwards import fwCS, fwMT, jointFW, tic, toc
from inverses_new import MTinv, CSinv,Jointinv
from Joint import Jointinv2, tar
import numpy.random as rnd
from scipy.spatial import distance
from scipy.optimize import minimize as smin
from scipy.signal import savgol_filter 

optmeth = 'Nelder-Mead'
m_iter = 300
regul = 0.1

def plotr(ro_and_h, style, col, lab):
    res=ro_and_h[0:int(len(ro_and_h)/2)+1]
    thick=ro_and_h[int(len(ro_and_h)/2)+1:]
    z=np.zeros((np.shape(res)))
    z[0]=0
    for i in range(0,len(thick)):
        z[i+1]=z[i]+thick[i]
    z=np.append(z,z[-1]+thick[-1])
    res=np.append(res[0],res)
    plt.step(z,res, style, color=col, label=lab)
    plt.legend()
    plt.yscale('log')
    
#_____________________________generate syntetic data___________________________
f = np.logspace(np.log10(10e-5), np.log10(10e+4), 100)
r = [30,25000,10,5000,2.5,1000]
h = [100,100,100,250,500]

mod=r+h

cs_field=(fwCS(f,mod))
mt_field =fwMT(f,mod,0)
obs_data=jointFW(f, mod)

rapr = [30,10000,10,6000,10,1000]
hapr = [50,150,100,300,500]

apr_mod = rapr+hapr
ref_mod = np.copy(apr_mod)



a = MTinv(apr_mod , f, mt_field, ref_mod, regul, m_iter, optmeth)
b = CSinv(apr_mod , f, cs_field, ref_mod, regul, m_iter, optmeth)
c = Jointinv(apr_mod , f, mt_field, cs_field, ref_mod, regul, m_iter, optmeth)
cc=Jointinv2(apr_mod , f, obs_data, ref_mod, regul, m_iter, optmeth)


plt.figure()
plt.loglog(1/f,fwMT(f,a.x,0),'--',color='red',label='MT')
plt.loglog(1/f,fwMT(f,c.x,0),'*',color='red',label='MT Joint2')
plt.loglog(1/f,fwMT(f,cc.x,0),'*',color='red',label='MT Joint2')
plt.loglog(1/f,mt_field,color='red')
plt.loglog(1/f,fwCS(f,b.x),'--',color='blue',label='CSEM')
plt.loglog(1/f,fwCS(f,c.x),'*',color='blue',label='CSEM Joint2')
plt.loglog(1/f,fwCS(f,cc.x),'*',color='blue',label='CSEM Joint2')
plt.loglog(1/f,cs_field,color='blue' )
plt.legend()
plt.xlabel('T, s')
plt.ylabel('App Resistivity, Ohm*m')
plt.show()

plt.figure()
plt.semilogy(obs_data)
plt.semilogy(jointFW(f, cc.x),'--',color='green')
plt.semilogy(jointFW(f, c.x),'--',color='purple')


plt.figure()
plt.xlabel('Depth, m')
plt.ylabel('Resistivity, Ohm*m')
plotr(mod,'-','black','True')
plotr(apr_mod,'--','pink','Apr')
plotr(a.x,'--','red','MT')
plotr(b.x,'--','blue','CSEM')
plotr(c.x,'--','purple','Joint1')
plotr(cc.x,'--','green','Joint2')
plt.show()














