# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 22:17:50 2020

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from forwards import fwCS, fwMT

def msf_mt(obj, data, f, ref, regul):
    curl=np.zeros(len(obj.allvecs))
    for i in range(len(obj.allvecs)):
        mod = obj.allvecs[i]
        syn_data=fwMT(f, mod, 0)    
            
        curF=np.sqrt(np.sum((((np.log(syn_data)-np.log(data))/np.log(data))**2)/len(data)))
        curM=np.sqrt(np.sum((((np.log(mod)-np.log(ref))/np.log(mod))**2)/len(mod)))
        
    #    curM=0
        curl[i] = curF + regul*curM
      
    return curl

def msf_cs(obj, data, f, ref, regul):
    curl=np.zeros(len(obj.allvecs))
    for i in range(len(obj.allvecs)):
        mod = obj.allvecs[i]
        syn_data=fwCS(f, mod)    
            
        curF=np.sqrt(np.sum((((np.log(syn_data)-np.log(data))/np.log(data))**2)/len(data)))
        curM=np.sqrt(np.sum((((np.log(mod)-np.log(ref))/np.log(mod))**2)/len(mod)))
        
    #    curM=0
        curl[i] = curF + regul*curM
    
    return curl

def msf_j(obj, mt_data, cs_data, f, ref, regul):
    curl=np.zeros(len(obj.allvecs))
    
    for i in range(len(obj.allvecs)):
        mod = obj.allvecs[i]
        
        synMT=fwMT(f, mod,0)
        curFMT = np.sqrt(np.sum((((np.log(synMT)-np.log(mt_data))/np.log(mt_data))**2)/len(mt_data)))
       
        synCS=fwCS(f, mod)
        curFCS = np.sqrt(np.sum((((np.log(synCS)-np.log(cs_data))/np.log(cs_data))**2)/len(cs_data)))
        
    
        curM = np.sqrt(np.sum((((np.log(mod)-np.log(ref))/np.log(ref))**2)/len(ref)))

        
    #    curM=0
        curl[i] = 0.5*curFMT+0.5*curFCS+regul*curM
    
    return curl