# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 21:17:48 2020

@author: User
"""

import numpy as np
from forwards import jointFW
from scipy.optimize import minimize as smin
Nfeval = 1

def callbackF_joint(Xi):
    global Nfeval
    print('Joint2 iter= ', Nfeval)
#    print('{0:4d}   {1: 3.2f}   {2: 3.2f}   {3: 3.2f}   {4: 3.2f}   {5: 3.2f}   {6: 3.2f}   {7: 3.2f}   {8: 3.2f}   {9: 3.2f}   {10: 3.2f}    {11: 3.2f}'.format(Nfeval2, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], Xi[5], Xi[6], Xi[7], Xi[8], Xi[9], Xi[10]))
    Nfeval += 1

def tar(cmod, f, real_data, ref, regul):
    syn_data=jointFW(f,cmod)
    curD=np.sum((np.log(np.array(syn_data))-np.log(np.array(real_data)))**2/np.log(np.array(real_data)))/len(real_data)
    curM=np.sum((np.log(np.array(cmod))-np.log(np.array(ref)))**2/np.log(np.array(ref)))/len(ref)    
    cur=np.sqrt(abs(curD))+regul*np.sqrt(abs(curM))
    
    for i in range(0, len(cmod)):
        if cmod[i]<=3:
            cur=999999
    
    return cur

def Jointinv2(apr_mod , f, real_data, ref_mod, regul, m_iter,optmeth):
    res = smin(tar,apr_mod,callback=callbackF_joint, args=(f, real_data, ref_mod, regul) , method = optmeth, options={'maxiter': m_iter} )
    return res