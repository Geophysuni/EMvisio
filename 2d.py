# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 18:21:21 2023

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio as iio

asis_image = iio.imread("oblaka-1.jpg")
ini_image = iio.imread("oblaka-1.jpg")[:, :, 0]


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
        
        xInd = int((tmpXCrd-tmpModX[0])/self.model.dx)
        iniModel = self.model.model[:,xInd]
        
        thickness = np.ones(int(self.model.shape[0]/self.ystep))*self.ystep
        resistivity = np.zeros(len(thickness)+1)
        
        for i in range(len(thickness)):
            resistivity[i] = np.mean(iniModel[int(i*self.ystep):int((i+1)*self.ystep)])
        resistivity[-1] = iniModel[-1]
    
        
        return thickness, resistivity
        

model = Model(ini_image)
model.setSize([0,1200], [0,338])

prb = FProblem('mt')
prb.setModel(model)
prb.setRec(np.arange(100,500,50))

# prb.setYStep(50)

# h, r = prb.get1DModel(1199)
