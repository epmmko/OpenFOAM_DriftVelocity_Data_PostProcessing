#### import the simple module from the paraview
from __future__ import division
from paraview import vtk
from paraview.simple import *
location=[1.6, 3.62, 6.81]
#Setting Up
lastTime=25
index=range(5*lastTime+1)
case=GetActiveSource()
#Slice
slice=Slice(Input=case)
slice.SliceType='Plane'
for l in location:
    time=[]
    j_a_time=[]
    a_d_time=[]
    v_d_time=[]
    #Slice
    slice.SliceType.Origin=[l, 0, 0.009525]
    slice.SliceType.Normal=[1, 0, 0]
    slice.Triangulatetheslice=0
    slice.Crinkleslice=0
    
    #MeshQuality
    meshq=MeshQuality(Input=slice)
    meshq.TriangleQualityMeasure='Area'
    meshq.QuadQualityMeasure='Area'
    
    for k in index:
       meshq.UpdatePipeline(k/5)
       data=paraview.servermanager.Fetch(meshq)
       numCells=data.GetBlock(0).GetBlock(0).GetNumberOfCells()
       cellData=data.GetBlock(0).GetBlock(0).GetCellData()
       #numCells=data.GetNumberOfCells()
       #cellData=data.GetCellData()
       
       alpha=cellData.GetArray('alpha.oil')
       U=cellData.GetArray('U')
       A=cellData.GetArray('Quality')
       
       Atotal=0
       a_d=0
       j_a=0
       j_d=0
       for i in range(numCells):
           Atotal=Atotal+A.GetValue(i)
           a_d=a_d+(1-alpha.GetValue(i))*A.GetValue(i)
           j_a=j_a+U.GetComponent(i,0)*A.GetValue(i)
           j_d=j_d+(1-alpha.GetValue(i))*U.GetComponent(i,0)*A.GetValue(i)
       a_d=a_d/Atotal
       j_a=j_a/Atotal
       j_d=j_d/Atotal
       time.append(k/5)
       a_d_time.append(a_d)
       j_a_time.append(j_a)
       if a_d==0:
           v_d_time.append(0)
       else:
           v_d_time.append(j_d/a_d)
       print('Location at '+str(l)+' and time is '+str(k/5))
#
    T=vtk.vtkTable()
    
    t=vtk.vtkDoubleArray()
    t.SetName("Time [s]")
    for i in time:
        t.InsertNextValue(i)
    T.AddColumn(t)
    
    alpha=vtk.vtkDoubleArray()
    alpha.SetName("Alpha_d")
    for i in a_d_time:
        alpha.InsertNextValue(i)
    T.AddColumn(alpha)
    
    v_j=vtk.vtkDoubleArray()
    v_j.SetName("V_j [m/s]")
    for i in j_a_time:
        v_j.InsertNextValue(i)
    T.AddColumn(v_j)
    
    v_d=vtk.vtkDoubleArray()
    v_d.SetName("V_d [m/s]")
    for i in v_d_time:
        v_d.InsertNextValue(i)
    T.AddColumn(v_d)
    
    Result=TrivialProducer()
    filter_aux=Result.GetClientSideObject()
    filter_aux.SetOutput(T)
    Result.UpdatePipeline()
    SaveData(str(l)+'.csv')
