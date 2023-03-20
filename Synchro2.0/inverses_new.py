# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:31:16 2020

@author: zhuravlev.sd
"""

import numpy as np
import matplotlib.pyplot as plt
from forwards import fwCS, fwMT, jointFW, tic, toc
import numpy.random as rnd
from scipy.spatial import distance
from scipy.optimize import minimize as smin
Nfeval = 1
Nfeval2 = 1
Nfeval3 = 1
Nfeval4 = 1


def tarMT(mod, f, real_data, ref_mod, regul):         
    syn_data=fwMT(f, mod, 0)    
    
#    curF = np.sum(np.sqrt(((np.log(syn_data)-np.log(real_data))**2)/np.log(real_data)))  
#    curM = np.sum(np.sqrt(((np.log(mod)-np.log(ref_mod))**2)/np.log(ref_mod)))
    
    curF=np.sqrt(np.sum((((np.log(syn_data)-np.log(real_data))/np.log(syn_data))**2)/len(real_data)))
    curM=np.sqrt(np.sum((((np.log(mod)-np.log(ref_mod))/np.log(mod))**2)/len(mod)))
    
    curF=np.sqrt(np.sum(((((syn_data)-(real_data))/(syn_data))**2)/len(real_data)))
    curM=np.sqrt(np.sum(((((mod)-(ref_mod))/(mod))**2)/len(mod)))
    
#    curM=0
    cur = curF + regul*curM
    
    apr_res=ref_mod[0:int(len(ref_mod)/2)+1]
    apr_thick=ref_mod[int(len(ref_mod)/2)+1:] 
    
    cur_res=mod[0:int(len(mod)/2)+1]
    cur_thick=mod[int(len(mod)/2)+1:]
        
    for i in range(len(cur_thick)):
        if abs(cur_thick[i]-apr_thick[i])/apr_thick[i]>0.1:
            cur=999999*cur

        
    for i in range(0, len(mod)):
        if mod[i]<=0:
            cur=999999
    

    return cur

def tarCS(mod, f, real_data, ref_mod, regul):         
    
    syn_data=(fwCS(f, mod))  
    
#    curF = np.sum(np.sqrt(((np.log(syn_data)-np.log(real_data))**2)/np.log(real_data)))  
#    curM = np.sum(np.sqrt(((np.log(mod)-np.log(ref_mod))**2)/np.log(ref_mod)))
    
#    curF=np.sqrt(np.sum((((np.log(syn_data)-np.log(real_data))/np.log(syn_data))**2)/len(real_data)))
#    curM=np.sqrt(np.sum((((np.log(mod)-np.log(ref_mod))/np.log(mod))**2)/len(mod)))
    
    curF=np.sqrt(np.sum(((((syn_data)-(real_data))/(syn_data))**2)/len(real_data)))
    curM=np.sqrt(np.sum(((((mod)-(ref_mod))/(mod))**2)/len(mod)))
    
#    curM=0
    
    cur = curF + regul*curM
    
    apr_res=ref_mod[0:int(len(ref_mod)/2)+1]
    apr_thick=ref_mod[int(len(ref_mod)/2)+1:] 
    
    cur_res=mod[0:int(len(mod)/2)+1]
    cur_thick=mod[int(len(mod)/2)+1:]
        
    for i in range(len(cur_thick)):
        if abs(cur_thick[i]-apr_thick[i])/apr_thick[i]>0.1:
            cur=999999*cur

        
    for i in range(0, len(mod)):
        if mod[i]<=0:
            cur=999999
    
    return cur

def callbackF_MT(Xi):
    global Nfeval
    print('MT iter= ', Nfeval)
#    print('{0:4d}   {1: 3.2f}   {2: 3.2f}   {3: 3.2f}   {4: 3.2f}   {5: 3.2f}   {6: 3.2f}   {7: 3.2f}   {8: 3.2f}   {9: 3.2f}   {10: 3.2f}    {11: 3.2f}'.format(Nfeval, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], Xi[5], Xi[6], Xi[7], Xi[8], Xi[9], Xi[10]))
    Nfeval += 1

def callbackF_CS(Xi):
    global Nfeval2
    print('CS iter= ', Nfeval2)
#    print('{0:4d}   {1: 3.2f}   {2: 3.2f}   {3: 3.2f}   {4: 3.2f}   {5: 3.2f}   {6: 3.2f}   {7: 3.2f}   {8: 3.2f}   {9: 3.2f}   {10: 3.2f}    {11: 3.2f}'.format(Nfeval2, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], Xi[5], Xi[6], Xi[7], Xi[8], Xi[9], Xi[10]))
    Nfeval2 += 1
    
def callbackF_J(Xi):
    global Nfeval3
#    print('Joint iter= ', Nfeval3)
#    print('Joint iter= ', Nfeval3)
#    print('{0:4d}   {1: 3.2f}   {2: 3.2f}   {3: 3.2f}   {4: 3.2f}   {5: 3.2f}   {6: 3.2f}   {7: 3.2f}   {8: 3.2f}   {9: 3.2f}   {10: 3.2f}    {11: 3.2f}'.format(Nfeval2, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], Xi[5], Xi[6], Xi[7], Xi[8], Xi[9], Xi[10]))
    Nfeval3 += 1
    
def callbackF_joint(Xi):
    global Nfeval4
    print('Joint2 iter= ', Nfeval4)
#    print('{0:4d}   {1: 3.2f}   {2: 3.2f}   {3: 3.2f}   {4: 3.2f}   {5: 3.2f}   {6: 3.2f}   {7: 3.2f}   {8: 3.2f}   {9: 3.2f}   {10: 3.2f}    {11: 3.2f}'.format(Nfeval2, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], Xi[5], Xi[6], Xi[7], Xi[8], Xi[9], Xi[10]))
    Nfeval4 += 1

def MTinv(apr_mod , f, MT_noise, ref_mod, regul, m_iter,optmeth):
    resMT = smin(tarMT, apr_mod , callback=callbackF_MT, args=(f, MT_noise, ref_mod, regul) , method = optmeth, options={'maxiter': m_iter, 'adaptive': True, 'return_all': True} )
    return resMT

def CSinv(apr_mod , f, CS_noise, ref_mod, regul, m_iter,optmeth):
    resCS = smin(tarCS, apr_mod , callback=callbackF_CS, args=(f, CS_noise, ref_mod, regul) , method = optmeth, options={'maxiter': m_iter, 'adaptive': True, 'return_all': True} )
    return resCS

def tarJoint(mod, f, MT, CS, ref_mod, regul):
    synMT=fwMT(f, mod,0)
#    curFMT = np.sqrt(np.sum((((np.log(synMT)-np.log(MT))/np.log(MT))**2)/len(MT)))
    curFMT = np.sqrt(np.sum(((((synMT)-(MT))/(MT))**2)/len(MT)))
    
    
    synCS=fwCS(f, mod)
#    curFCS = np.sqrt(np.sum((((np.log(synCS)-np.log(CS))/np.log(CS))**2)/len(CS)))
    curFCS = np.sqrt(np.sum(((((synCS)-(CS))/(CS))**2)/len(CS)))
    

#    curM = np.sqrt(np.sum((((np.log(mod)-np.log(ref_mod))/np.log(ref_mod))**2)/len(ref_mod)))
    curM = np.sqrt(np.sum(((((mod)-(ref_mod))/(ref_mod))**2)/len(ref_mod)))
#    syndata=np.concatenate((synMT,synCS),axis=0)
#    data=np.concatenate((MT,CS),axis=0)
#    
#    curdat=np.sum(np.sqrt(((np.log(syndata)-np.log(data))**2)/np.log(data)))  
     
#    curM=0    
    w=curFMT/curFCS    #считаем, что curFMT<curFCS
    if w>=1:
        w=1/w        
    ww=1-w
    
    if curFMT>=curFCS:
        w1=np.max([w,ww])
        w2=np.min([w,ww])
    else:
        w1=np.min([w,ww])
        w2=np.max([w,ww])
#    w1,w2=1,1
        
    print('w1='+str(w1)+' '+'w2='+str(w2))
    
#    regul=0
    
    cur=w1*curFMT+w2*curFCS+regul*curM
    
#    cur=curdat+regul*curM
    
    apr_res=ref_mod[0:int(len(ref_mod)/2)+1]
    apr_thick=ref_mod[int(len(ref_mod)/2)+1:] 
    
    cur_res=mod[0:int(len(mod)/2)+1]
    cur_thick=mod[int(len(mod)/2)+1:]
        
    for i in range(len(cur_thick)):
        if abs(cur_thick[i]-apr_thick[i])/apr_thick[i]>0.1:
            cur=999999*cur
    
    for i in range(0, len(mod)):
        if mod[i]<=0:
            cur=999999
    
    return cur

def Jointinv(apr_mod , f, MT_noise, CS_noise, ref_mod, regul, m_iter,optmeth):        
    res = smin(tarJoint, apr_mod , callback=callbackF_J, args=(f, MT_noise, CS_noise, ref_mod, regul) , method = optmeth, options={'maxiter': m_iter, 'adaptive': True, 'return_all': True} )
    return res
#
#def tarJoint2(mod, f, obs_data, ref_mod, regul):
#    print(mod)
#    syn_data=jointFW(f, mod)
#    cur_data = np.sum(np.sqrt((np.log(syn_data)-np.log(obs_data))**2)/np.log(obs_data))/len(obs_data)
#        
#    curM = np.sum(np.sqrt((np.log(mod)-np.log(ref_mod))**2)/np.log(ref_mod))/len(ref_mod)
#    
#    regul=0
#    
#    cur=cur_data+regul*curM
#    
#    apr_res=ref_mod[0:int(len(ref_mod)/2)+1]
#    apr_thick=ref_mod[int(len(ref_mod)/2)+1:] 
#    
#    cur_res=mod[0:int(len(mod)/2)+1]
#    cur_thick=mod[int(len(mod)/2)+1:]
#        
#    for i in range(len(cur_thick)):
#        if abs(cur_thick[i]-apr_thick[i])/apr_thick[i]>0.1:
#            cur=999999*cur
#    
#    for i in range(0, len(mod)):
#        if mod[i]<=0:
#            cur=999999    
#    return cur
#
#def Jointinv2(apr_mod , f, obs_data, ref_mod, regul, m_iter,optmeth):
#    res = smin(tarJoint2, apr_mod, callback=callbackF_joint, args=(f, obs_data, ref_mod, regul) , method = optmeth, options={'maxiter': m_iter} )
#    return res































