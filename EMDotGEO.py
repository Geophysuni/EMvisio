# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 22:55:03 2023

@author: Sergey Zhuravlev
"""

from tkinter import *
from  tkinter import ttk
from tkinter import filedialog as fd

import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import pickle
import pandas as pd
import numpy as np
import imageio as iio
from scipy import interpolate
from emlibdotgeo import Model, FProblem
from importGRD import importGRD


mainWin = Tk()
mainWin.geometry('1300x720')
mainWin.title('FWDEM DotGEO')

# elements
createModelBut = Button(master=mainWin, width=25, text = 'Create model')
setObsSystemBut = Button(master=mainWin, width=25, text = 'Obs. system')

prbType = IntVar()
csemRB = Radiobutton(master=mainWin, text="CSEM", variable=prbType, value=1)
mtRB = Radiobutton(master=mainWin, text="MT", variable=prbType, value=2)

frqMinLabel = Label(master=mainWin, text = 'Min f:')
frqMinEntry = Entry(master=mainWin, width = 10)
frqMaxLabel = Label(master=mainWin, text = 'Max f:')
frqMaxEntry = Entry(master=mainWin, width = 10)
frqNLabel = Label(master=mainWin, text = 'n:')
frqNEntry = Entry(master=mainWin, width = 10)

calculateBut = Button(master=mainWin, width=25, text = 'Calculate')
exportBut = Button(master=mainWin, width=25, text = 'Export')

pointListBox = Listbox(master=mainWin, height = 28, width = 30)

fig_model = Figure(figsize = (10, 3), dpi = 100)
canvas_model = FigureCanvasTkAgg(fig_model, master = mainWin)

fig_measurements = Figure(figsize = (10, 3.5), dpi = 100)
canvas_measurements = FigureCanvasTkAgg(fig_measurements, master = mainWin)




def openCreateModel(event):
    createModLayer = Toplevel(mainWin)
    createModLayer.geometry("1100x600")
    
    # elements
    matFrame = LabelFrame(master = createModLayer, text = 'Matrix', height = 170, width = 210)
    
    importFromFileBut = Button(master=matFrame, width=25, text = 'Import from file')
    importFromLayer = Button(master=matFrame, width=25, text = 'Add layer')
    
    minxLabel = Label(master = matFrame, text = 'Min X')
    minxEntry = Entry(master = matFrame, width=5)
    maxxLabel = Label(master = matFrame, text = 'Max X')
    maxxEntry = Entry(master = matFrame, width=5)
    
    minzLabel = Label(master = matFrame, text = 'Min Z')
    minzEntry = Entry(master = matFrame, width=5)
    maxzLabel = Label(master = matFrame, text = 'Max Z')
    maxzEntry = Entry(master = matFrame, width=5)
    
    setSizeBut = Button(master=matFrame, width=25, text = 'Set sizes')
    
    saveBut = Button(master=createModLayer, width=25, text = 'Save')
    
    fig_model_win = Figure(figsize = (8, 5.1), dpi = 100)
    canvas_model_win = FigureCanvasTkAgg(fig_model_win, master = createModLayer)
    
    # functions
    def importFromMatr(event):
        filename = fd.askopenfilename(initialdir = os.getcwd())
        
        xmin, xmax, zmin, zmax, matrix = importGRD(filename)
        
        minxEntry.delete(0,END)
        maxxEntry.delete(0,END)
        minzEntry.delete(0,END)
        maxzEntry.delete(0,END)
        
        minxEntry.insert(END, str(xmin))
        maxxEntry.insert(END, str(xmax))
        minzEntry.insert(END, str(zmin))
        maxzEntry.insert(END, str(zmax))
        
        
        fig_model_win.clf()
        plt1 = fig_model_win.add_subplot(111)     
        
        img = plt1.imshow(matrix)
        fig_model_win.colorbar(img, label = 'R, Ohm*m')
        
        fig_model_win.tight_layout()
        canvas_model_win.draw()
        
        try:
            with open('tmpmd', 'rb') as f:
                modpar = pickle.load(f)
            modpar['matrix'] = matrix
            with open('tmpmd', 'wb') as f:
                pickle.dump(modpar, f)
        except:
            modpar = {}
            modpar['matrix'] = matrix
            with open('tmpmd', 'wb') as f:
                pickle.dump(modpar, f)
        
        return
    
    
    
    def setSize(event):
        
        xmin, xmax = np.float(minxEntry.get()), np.float(maxxEntry.get())
        zmin, zmax = np.float(minzEntry.get()), np.float(maxzEntry.get())
        
        with open('tmpmd', 'rb') as f:
            modpar = pickle.load(f)
            
        modpar['xlim'] = [xmin, xmax]
        modpar['zlim'] = [zmin, zmax]
        
        extent = (xmin, xmax, zmax, zmin)
        fig_model_win.clf()
        plt1 = fig_model_win.add_subplot(111)     
        
        img = plt1.imshow(modpar['matrix'], extent = extent, aspect='auto')
        fig_model_win.colorbar(img, label = 'R, Ohm*m')
        
        fig_model_win.tight_layout()
        canvas_model_win.draw()
        
        with open('tmpmd', 'wb') as f:
            pickle.dump(modpar, f)
        
    def saveMod(event):
        # with open('tmpmd', 'wb') as f:
        #     pickle.dump(modpar, f)
            
        with open('tmpmd', 'rb') as f:
            modpar = pickle.load(f)
        
        extent = (modpar['xlim'][0], modpar['xlim'][1], modpar['zlim'][1], modpar['zlim'][0])
        fig_model.clf()
        plt1 = fig_model.add_subplot(111)     
        
        img = plt1.imshow(modpar['matrix'], extent = extent, aspect='auto')
        fig_model.colorbar(img, label = 'R, Ohm*m')
        
        fig_model.tight_layout()
        canvas_model.draw()
        
        createModLayer.destroy()
        
        return
    # binding
    importFromFileBut.bind('<ButtonRelease-1>', importFromMatr)
    setSizeBut.bind('<ButtonRelease-1>', setSize)
    saveBut.bind('<ButtonRelease-1>', saveMod)
    # placing
    
    matFrame.place(x=10, y=10)
    importFromFileBut.place(x=10, y=10)
        
    minxLabel.place(x=10, y=40)
    minxEntry.place(x=50, y=40)
    maxxLabel.place(x=100, y=40)
    maxxEntry.place(x=140, y=40)
    
    minzLabel.place(x=10, y=70)
    minzEntry.place(x=50, y=70)
    maxzLabel.place(x=100, y=70)
    maxzEntry.place(x=140, y=70)
    
    setSizeBut.place(x=10, y=100)
    
    listFrame = LabelFrame(master = createModLayer, text = 'layers', height = 370, width = 210)
    
    layListFrame = ttk.Treeview(master = listFrame, height=10)
    layListFrame['columns'] = ('Layer','Res')
    layListFrame.column('#0', width=0, stretch=NO)
    layListFrame.column('Layer', width=90)
    layListFrame.column('Res', width=90)
    layListFrame.heading('Layer', text='Layer')
    layListFrame.heading('Res', text='Res') 
    
    importLayerFileBut = Button(master=listFrame, width=25, text = 'Import from file')
    delLayerBut = Button(master=listFrame, width=25, text = 'Delete')
    setLayerBut = Button(master=listFrame, width=15, text = 'Set value')
    rValueEntry = Entry(master=listFrame, width=8)
    
    
    listFrame.place(x=10, y=160)
    importLayerFileBut.place(x=10, y=10)
    delLayerBut.place(x=10, y=40)
    setLayerBut.place(x=10, y=70)
    rValueEntry.place(x=130, y=72)
    layListFrame.place(x=10, y=100)
    
    
    def plotLays():
        
        
        
        try:
                        
            layMatr = []
            
            df = pd.read_csv('laylist.csv', sep = '/', index_col=0)
            for i, path in enumerate(list(df['path'])):
                with open(path, 'r') as f:
                    clay = f.readlines()
                cparr = []
                for cp in clay:
                    cparr.append(cp.replace('\n', '').split('	'))
                layMatr.append((np.array(cparr).astype(float), list(df['val'])[i]))
                
            
            minX, maxX = 99999, -99999
            for el in layMatr:
                if np.min(el[0][:,0])<minX:
                    minX = np.min(el[0][:,0])
                if np.max(el[0][:,0])>maxX:
                    maxX = np.max(el[0][:,0])
            
            newX = np.arange(minX, maxX)
            newLays = []
            for el in layMatr:
                f = interpolate.interp1d(el[0][:,0], el[0][:,1], fill_value='extrapolate')
                newY = f(newX)
                newLays.append((np.max(newY), newY, el[1]))
            newLays.sort()
            
            model = np.zeros((int(newLays[-1][0]), len(newX)))
            for cl in newLays:
                for j in range(len(newX)):
                    model[int(cl[1][j]):, j] = cl[2]
                    
            
            fig_model_win.clf()
            plt1 = fig_model_win.add_subplot(111)
            extent = [minX, maxX, newLays[-1][0], np.min(newLays[0][1])]
            img = plt1.imshow(model, extent = extent, aspect='auto')
            for lay in newLays:
                plt1.plot(newX, lay[1], color = 'black')
            fig_model_win.colorbar(img, label = 'R, Ohm*m')
           
            fig_model_win.tight_layout()
            canvas_model_win.draw()
            
            modpar = {}
            modpar['xlim'] = [minX, maxX]
            modpar['zlim'] = [np.min(newLays[0][1]), newLays[-1][0]]
            modpar['matrix'] = model
            
            with open('tmpmd', 'wb') as f:
                pickle.dump(modpar, f)
                    
                
        except:
            bonk = 1
    
    def importLayer(event):
        


            
        if len(layListFrame.get_children())==0:
            
            filenames = fd.askopenfilenames(initialdir = os.getcwd())
                    
            i=0
            for file in filenames:
                fname = file.split('/')[-1].split('.')[0]
                layListFrame.insert('', i, values=(fname))
                print(fname)
                i = i+1
            
            df = pd.DataFrame(columns = ['path', 'name', 'val'])
            for i, file in enumerate(filenames):
                fname = file.split('/')[-1].split('.')[0]
                df.loc[i,'path'] = file
                df.loc[i,'name'] = fname
        else:
            filenames = fd.askopenfilenames(initialdir = os.getcwd())
                    
            i=0
            for file in filenames:
                fname = file.split('/')[-1].split('.')[0]
                layListFrame.insert('', i, values=(fname))
                print(fname)
                i = i+1
                
            df = pd.read_csv('laylist.csv', sep = '/', index_col=0).reset_index(drop=True)
            iniLen = len(df)
            i = 0 
            for file in filenames:
                fname = file.split('/')[-1].split('.')[0]
                if fname not in list(df['name']):
                    df.loc[iniLen+i,'path'] = file
                    df.loc[iniLen+i,'name'] = fname
                    i=i+1
        
        
        # if layListFrame.size()>0:
            
        #     print('importLayer initial')
        #     df = pd.read_csv('laylist.csv', sep = '/', index_col=0).reset_index(drop=True)
        #     iniLen = len(df)
        #     i = 0 
        #     for file in filenames:
        #         fname = file.split('/')[-1].split('.')[0]
        #         if fname not in list(df['name']):
        #             df.loc[iniLen+i,'path'] = file
        #             df.loc[iniLen+i,'name'] = fname
        #             i=i+1
            
        # else:
        #     print('importLayer exception')
        #     df = pd.DataFrame(columns = ['path', 'name', 'val'])
        #     for i, file in enumerate(filenames):
        #         fname = file.split('/')[-1].split('.')[0]
        #         df.loc[i,'path'] = file
        #         df.loc[i,'name'] = fname
        
        df.to_csv('laylist.csv', sep = '/')
        
        plotLays()
        
    def setValue(event):
        
        curVal = np.float(rValueEntry.get())
        
        selection = layListFrame.focus()
        index = layListFrame.index(selection)
        
        name, val = layListFrame.item(selection)['values'][0], curVal
        layListFrame.insert('', index, values=(name, val))
        layListFrame.delete(selection)
        
        df = pd.read_csv('laylist.csv', sep = '/', index_col=0).reset_index(drop=True)
        
        
        ind = df[df['name']==name].index[0]
        df.loc[ind, 'val'] = val
        df.to_csv('laylist.csv', sep = '/')
        
        plotLays()  
        
    def delLayer(event):
        
        selected_item = layListFrame.selection()[0] 
        layListFrame.delete(selected_item)
        
        updatedList = []
        treeList = layListFrame.get_children()
        restNames = [] 
        for el in treeList:
            restNames.append(layListFrame.item(el)["values"][0])
        df = pd.read_csv('laylist.csv', sep = '/', index_col=0)
        newDf = df[df['name'].isin(restNames)]
        newDf.to_csv('laylist.csv', sep = '/')
        plotLays()
            
        
        
        
    
    importLayerFileBut.bind('<ButtonRelease-1>', importLayer)
    setLayerBut.bind('<ButtonRelease-1>', setValue)
    delLayerBut.bind('<ButtonRelease-1>', delLayer)
    
    saveBut.place(x=20, y=540)
    
    canvas_model_win.get_tk_widget().place(x = 250, y = 20)
    
    createModLayer.mainloop()

# observationsystem
def openCreateObs(event):
    createObsLayer = Toplevel(mainWin)
    createObsLayer.geometry("1200x550")
    
    
    importBut = Button(master=createObsLayer, width=25, text = 'Import profile')
    addPointBut = Button(master=createObsLayer, width=15, text = 'Add point')
    pointEntry = Entry(master=createObsLayer, width=10)
    delPointBut = Button(master=createObsLayer, width=25, text = 'Delete point')
    saveBut = Button(master=createObsLayer, width=25, text = 'Save')
    
    fig_model_obs = Figure(figsize = (9, 5), dpi = 100)
    canvas_model_obs = FigureCanvasTkAgg(fig_model_obs, master = createObsLayer)
    
    pointList = Listbox(master=createObsLayer, width = 30, height = 23)
    
    def plotObs():
        with open('points', 'rb') as f:
            points = pickle.load(f)
        
        with open('tmpmd', 'rb') as f:
            modpar = pickle.load(f)
         
        fig_model_obs.clf()
        plt1 = fig_model_obs.add_subplot(111)
        
        extent = [modpar['xlim'][0], modpar['xlim'][1], modpar['zlim'][1], modpar['zlim'][0]]
        img = plt1.imshow(modpar['matrix'], extent = extent, aspect='auto')
        fig_model_obs.colorbar(img, label = 'R, Ohm*m')
        
        try:
            plt1.axvline(x=points, color = 'red')
        except:
            for point in points:
                plt1.axvline(x=point, color = 'red')
        
        fig_model_obs.tight_layout()
        canvas_model_obs.draw()
        
    
    
    def importProfile(event):
        
        filename = fd.askopenfilename(initialdir = os.getcwd())
        
        with open(filename, 'r') as f:
            points = np.array(f.readlines()).astype(float)
        for point in points:
            pointList.insert(END, point) 
        with open('points', 'wb') as f:
            pickle.dump(points, f)
        plotObs()
    
    def addPoint(event):
        
        allPoints = pointList.get(0, END)
        if len(allPoints)==0:
            newP = float(pointEntry.get())
            points = np.array(newP)
            pointList.insert(END, newP)
        else:
            try:
                with open('points', 'rb') as f:
                    points = pickle.load(f)
                newP = float(pointEntry.get())
                points = np.hstack((points, newP))
                pointList.insert(END, newP)
            except:
                points = np.array(newP)
        with open('points', 'wb') as f:
            pickle.dump(points, f)
        plotObs()
    
    def delPoint(event):
        selection = pointList.curselection()        
        pointList.delete(selection)
        
        allPoints = pointList.get(0, END)
        with open('points', 'wb') as f:
            pickle.dump(allPoints, f)
          
        plotObs()
        
    def saveObs(event):
        
        with open('points', 'rb') as f:
            points = pickle.load(f)
            
        with open('tmpmd', 'rb') as f:
            modpar = pickle.load(f)
        
        extent = (modpar['xlim'][0], modpar['xlim'][1], modpar['zlim'][1], modpar['zlim'][0])
        fig_model.clf()
        plt1 = fig_model.add_subplot(111)     
        
        img = plt1.imshow(modpar['matrix'], extent = extent, aspect='auto')
        fig_model.colorbar(img, label = 'R, Ohm*m')
        
        for point in points:
            plt1.axvline(x=point, color = 'red')
        
        fig_model.tight_layout()
        canvas_model.draw()
        
        pointListBox.delete(0,END)
        for point in points:
            pointListBox.insert(END,point)
        
        createObsLayer.destroy()
    
    importBut.bind('<ButtonRelease-1>', importProfile)
    addPointBut.bind('<ButtonRelease-1>', addPoint)
    delPointBut.bind('<ButtonRelease-1>', delPoint)
    saveBut.bind('<ButtonRelease-1>', saveObs)
    
    importBut.place(x=10,y=10)
    addPointBut.place(x=10,y=40)
    pointEntry.place(x=130,y=42)
    delPointBut.place(x=10,y=70)
    pointList.place(x=10,y=100)
    saveBut.place(x=10,y=480)
    canvas_model_obs.get_tk_widget().place(x = 250, y = 10)
    
    createObsLayer.mainloop()
    
 
def runCalc(event):
    with open('tmpmd', 'rb') as f:
        modpar = pickle.load(f)
    with open('points', 'rb') as f:
        obsPoints = pickle.load(f)
    
    if prbType.get() == 1:
        prob = 'cs'
    else:
        prob = 'mt'
    
    f = np.logspace(np.log10(float(frqMinEntry.get())), np.log10(float(frqMaxEntry.get())), int(frqNEntry.get()))
    
    f = np.logspace(np.log10(5*10e-5), np.log10(10e+3), 120)
    
    model = Model(modpar['matrix'])
    model.setSize(modpar['xlim'], modpar['zlim'])

    prb = FProblem(prob)
    prb.setModel(model)
    prb.setMeasurRenge(f)
    prb.setRec(obsPoints)
    
    prb.calculate()
    
    with open('curPrb','wb') as f:
        pickle.dump(prb, f)

def plotMeasure(event):
    
    selection = pointListBox.curselection()[0]
    
    
    with open('curPrb','rb') as f:
        prb = pickle.load(f)
    
    fig_measurements.clf()
    plt1 = fig_measurements.add_subplot(111)     
    
    
    plt1.loglog(prb.measureRange, prb.measuredValues[selection,:], color = 'black')
    plt1.set_xlabel('frequency')
    plt1.set_ylabel('App.resistivity, Ohm*m')
    
    fig_measurements.tight_layout()
    canvas_measurements.draw()
    
    with open('tmpmd', 'rb') as f:
        modpar = pickle.load(f)
    with open('points', 'rb') as f:
        points = pickle.load(f)
        
    extent = (modpar['xlim'][0], modpar['xlim'][1], modpar['zlim'][1], modpar['zlim'][0])
    fig_model.clf()
    plt1 = fig_model.add_subplot(111)     
    
    img = plt1.imshow(modpar['matrix'], extent = extent, aspect='auto')
    fig_model.colorbar(img, label = 'R, Ohm*m')
    
    
    for ip, point in enumerate(points):
        width = 1
        if ip==selection:
            width = 3
        plt1.axvline(x=point, color = 'red', linewidth = width)
    
    fig_model.tight_layout()
    canvas_model.draw()
         

def exportMeasurements(event):
    
    with open('curPrb','rb') as f:
        prb = pickle.load(f)
        
    columns = ['freq']
    for el in prb.rec:
        columns.append('point_'+str(el))
        
    file = fd.asksaveasfilename(defaultextension=".csv",filetypes=[("csv file", ".csv")])
    expDF = pd.DataFrame(columns = columns)
    expDF['freq'] = prb.measureRange
    icol = 1
    for el in prb.measuredValues:
        expDF[columns[icol]] = el
        icol = icol +1
    expDF.to_csv(file, sep = '\t')   
    

# binding
createModelBut.bind('<ButtonRelease-1>', openCreateModel)
setObsSystemBut.bind('<ButtonRelease-1>', openCreateObs)
calculateBut.bind('<ButtonRelease-1>', runCalc)
pointListBox.bind('<ButtonRelease-1>', plotMeasure)
exportBut.bind('<ButtonRelease-1>', exportMeasurements)

# placing
createModelBut.place(x=10,y=10)
setObsSystemBut.place(x=10,y=40)
csemRB.place(x=30,y=70)
mtRB.place(x=100,y=70)

frqMinLabel.place(x=10,y=100)
frqMinEntry.place(x=50,y=100)
frqMaxLabel.place(x=10,y=130)
frqMaxEntry.place(x=50,y=130)
frqNLabel.place(x=10,y=160)
frqNEntry.place(x=50,y=160)

calculateBut.place(x = 10, y = 190)
exportBut.place(x = 10, y = 220)

pointListBox.place(x = 10, y = 250)

canvas_measurements.get_tk_widget().place(x = 250, y = 10)
canvas_model.get_tk_widget().place(x = 250, y = 400)

# mainWin.after(20, update)
mainWin.mainloop()