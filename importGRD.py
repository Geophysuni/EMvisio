# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 18:52:21 2023

@author: Sergey Zhuravlev
"""

import numpy as np
# import matplotlib.pyplot as plt

def importGRD(file):

    with open(file, 'r') as f:
        inidat = f.readlines()
    
    tmpData = ''
    for string in inidat:
        tmpData = tmpData+string.replace('\n','')
    
    data = []
    tmp = ''
    for i, el in enumerate(tmpData):
        if el!=' ':
            tmp = tmp+el
        else:
            if len(tmp)!=0:
                data.append(tmp)
                tmp = ''
                
    nx, nz = int(data[1]),int(data[2])
    xmin, xmax = float(data[3]),float(data[4])
    zmin, zmax = float(data[5]),float(data[6]) 
    # marnge = [float(data[7]),float(data[8])]
    matrix = np.array(data[9:]).astype(float)
    
    matrix = matrix.reshape(nz,nx)

    return xmin, xmax, zmin, zmax, matrix+100
