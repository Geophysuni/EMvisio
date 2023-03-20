# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 19:17:09 2020

@author: zhuravlev.sd
"""


import numpy as np
import matplotlib.pyplot as plt
from forwards import fwCS, fwMT, jointFW, tic, toc
from inverses_new import MTinv, CSinv, Jointinv #, Jointinv2
import numpy.random as rnd
from scipy.spatial import distance
from scipy.optimize import minimize as smin
from scipy.signal import savgol_filter

from misfits import msf_mt, msf_cs, msf_j 

optmeth = 'Nelder-Mead'
m_iter = 250
regul = 0.2

#__________________________________plot stuff__________________________________

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

def mis_matrix(m,m1,m2,m3):
    r=m[0:int(len(m)/2)+1]
    h=m[int(len(m)/2)+1:]
    
    r1=m1[0:int(len(m1)/2)+1]
    h1=m1[int(len(m1)/2)+1:]
    
    r2=m2[0:int(len(m2)/2)+1]
    h2=m2[int(len(m2)/2)+1:]
    
    rj=m3[0:int(len(m3)/2)+1]
    hj=m3[int(len(m3)/2)+1:]
    
    matrix=np.zeros((len(r),14))
    matrix[:,0]=r
    matrix[:-1,1]=h
    
    matrix[:,2]=r1
    matrix[:-1,3]=h1
    
    matrix[:,4]=r2
    matrix[:-1,5]=h2
    
    matrix[:,6]=rj
    matrix[:-1,7]=hj
    
    #mist in %
    
    matrix[:,8]=100*(np.array(r)-np.array(r1))/np.array(r)
    matrix[:-1,9]=100*(np.array(h)-np.array(h1))/np.array(h)
    matrix[:,10]=100*(np.array(r)-np.array(r2))/np.array(r)
    matrix[:-1,11]=100*(np.array(h)-np.array(h2))/np.array(h)
    matrix[:,12]=100*(np.array(r)-np.array(rj))/np.array(r)
    matrix[:-1,13]=100*(np.array(h)-np.array(hj))/np.array(h)
    
    dev=np.zeros((6))
    
    dev[0]=np.sqrt(np.sum(matrix[:,8]/len(r))**2)
    dev[1]=np.sqrt(np.sum(matrix[:-1,9]/len(h))**2)
    
    dev[2]=np.sqrt(np.sum(matrix[:,10]/len(r))**2)
    dev[3]=np.sqrt(np.sum(matrix[:-1,11]/len(h))**2)
    
    dev[4]=np.sqrt(np.sum(matrix[:,12]/len(r))**2)
    dev[5]=np.sqrt(np.sum(matrix[:-1,13]/len(h))**2)
    
    
#    print(r)
    
    
    return matrix, dev


#_____________________________generate syntetic data___________________________
f = np.logspace(np.log10(5*10e-5), np.log10(10e+3), 120)
r = [30,25000,10,5000,2.5,1000]
h = [100,100,100,250,500]

#r=[20,9,15,30,1,60,10,1000]
#h=[100,260,240,230,160,270,1000]
mod=r+h

CSDATA=fwCS(f,mod)
MTDATA_r =fwMT(f,mod,0)
#_____________________________ad noise_________________________________________
noise=0.100
MT_noise=np.zeros((np.shape(MTDATA_r)))
CS_noise=np.zeros((np.shape(CSDATA)))

for i in range(0,len(MTDATA_r)):
    tmp=rnd.random()
    MT_noise[i]=MTDATA_r[i]+MTDATA_r[i]*(tmp-0.5)*2*noise


for i in range(0,len(CS_noise)):
    tmp=rnd.random()
    CS_noise[i]=CSDATA[i]+CSDATA[i]*(tmp-0.5)*2*noise




#____________________________smooth noise______________________________________

mt_field=savgol_filter(MT_noise, 11, 5)
cs_field=savgol_filter(CS_noise, 11, 5)


#mt_field=MT_noise
#cs_field=CS_noise

#obs_data=np.append(mt_field,cs_field)

#____________________________go inverse________________________________________

#r = [30,25000,10,5000,2.5,1000]
#h = [100,100,100,250,500]

#r = [30,25000,10,5000,2.5,1000]
#h = [100,100,100,250,500]

rapr = [30,10000,5,7000,10,1200]
hapr = [100,100,100,250,500]

#rapr=[30,10,20,30,2,50,5,1200]
#hapr=[100,260,240,230,160,270,1000]

apr_mod = rapr+hapr
ref_mod = rapr+hapr


#joint_data=jointFW(f,ref_mod)

#
a = MTinv(apr_mod , f, mt_field, ref_mod, regul, m_iter, optmeth)
b = CSinv(apr_mod , f, cs_field, ref_mod, regul, m_iter, optmeth)
c = Jointinv(apr_mod , f, mt_field, cs_field, ref_mod, regul, m_iter, optmeth)
#c=Jointinv2(apr_mod , f, joint_data, ref_mod, regul, m_iter, optmeth)
#a1,b1 = Jointinv(apr_mod , f, mt_field, cs_field, ref_mod, regul, m_iter, optmeth)



matrix,dev = mis_matrix(mod,a.x,b.x,c.x)

plt.figure()
plt.loglog(1/f,fwMT(f,a.x,0),'--',color='red',label='МТЗ')
plt.loglog(1/f,mt_field,color='red')
plt.loglog(1/f,fwCS(f,b.x),'--',color='blue',label='ЧЗ')
plt.loglog(1/f,cs_field,color='blue' )
plt.legend()
plt.xlabel('T, сек')
plt.ylabel('Каж. сопротивление, Ом*м')
plt.show()

plt.figure()
plt.loglog(1/f,fwCS(f,c.x),'--',color='blue',label='ЧЗ')
plt.loglog(1/f,cs_field,color='blue' )
plt.loglog(1/f,fwMT(f,c.x,0),'--',color='red',label='МТЗ')
plt.loglog(1/f,mt_field,color='red')
plt.legend()
plt.xlabel('T, сек')
plt.ylabel('Каж. сопротивление, Ом*м')
plt.show()

plt.figure()
plt.xlabel('Глубина, м')
plt.ylabel('УЭС, Ом*м')
plotr(mod,'-','black','Ист.')
plotr(apr_mod,'--','grey','Старт.')
plotr(a.x,'--','red','по МТЗ')
plotr(b.x,'--','blue','по ЧЗ')
plotr(c.x,'--','green','Синхр.')
plt.show()

