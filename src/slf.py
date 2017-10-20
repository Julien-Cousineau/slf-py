#! /usr/bin/python
import os,sys
import numpy as np
from fromOpenTelemac.parserSELAFIN import SELAFIN

class SLF(SELAFIN):
  def __init__(self,filePath=''):
    self.empty = True
    self.filePath = filePath
    SELAFIN.__init__(self, filePath)
    self._extent = None
    self._trixy = None
    self._triarea = None
    self._tricentroid = None
    self._values = None        
  
  @property
  def filePath(self):
    return self._filePath

  @filePath.setter
  def filePath(self, x):
    if(x != ''):
      if not os.path.exists(x): sys.exit("WARNING : Selafin path ({0}) does not exist".format(x))
      self.empty = False
    self._filePath = x

  @property
  def values(self):
    if self._values is None:
      values = np.zeros((self.NFRAME,self.NVAR,self.NPOIN3),np.float32)
      if not self.empty:
        for t in range(self.NFRAME):
          values[t] = self.getVALUES(t)
      self._values = values
    return self._values
      
  @values.setter
  def values(self,x):
    if x.shape[1] != self.NVAR: sys.exit("WARNING : Please check # of variables in setValues")
    self._values = x
    self.tags['times'] = np.arange(x.shape[0])
    
  @property
  def EXTENT(self):
    if self._extent is None:
      self._extent = np.array([np.min(self.MESHX),np.min(self.MESHY),np.max(self.MESHX),np.max(self.MESHY)])
    return self._extent
    

  @property
  def NFRAME(self):
    return len(self.tags['times'])
  
  @property
  def trixy(self):
    if self._trixy is None:
        x,y=(self.MESHX,self.MESHY)
        z = np.zeros(x.size)
        indexes = self.IKLE2
        self._trixy = np.concatenate((x[indexes][...,np.newaxis],y[indexes][...,np.newaxis],z[indexes][...,np.newaxis]), axis=2)
    return self._trixy
  
  @property
  def triarea(self):
    if self._triarea is None:
      trixy = self.trixy
      v0 = trixy[:,1,] - trixy[:,0,:]
      v1 = trixy[:,2,:] - trixy[:,0,:]
      self._triarea = 0.5* np.cross(v0,v1)[:,2]
    return self._triarea

  @property
  def tricentroid(self):
    if self._tricentroid is None:
      self._tricentroid = np.mean(self.trixy,axis=1)
    return self._tricentroid

   # {STRING} title
  def addTITLE(self,title):
    self.TITLE = '{: <{}}'.format(title, 80)
   
   # {OBJECT (name:str,unit:str)} var
  def addVAR(self,var):
    self.NBV1 += 1
    self.NVAR = self.NBV1 + self.NBV2
    self.VARINDEX = range(self.NVAR)
    self.VARNAMES.append('{: <{}}'.format(var['name'], 16)); 
    self.VARUNITS.append('{: <{}}'.format(var['unit'], 16));
    
  # {2D Array}
  def addPOIN(self,points):
    self.IPOB2 = np.arange(len(points))
    self.IPOB3 = self.IPOB2
    self.IPARAM = np.zeros(10)
    self.IPARAM[0] = 1
    self.NPOIN2 = len(points)
    self.NPOIN3 =self.NPOIN2
    self.MESHX = points[:,0];
    self.MESHY = points[:,1];

  # {2D Array(NELEM,3}
  def addIKLE(self,ikle):
    self.NDP2 = 3
    self.NDP3 = 3
    self.NELEM2 = len(ikle)
    self.NELEM3 = self.NELEM2
    self.IKLE2 = ikle
    self.IKLE3 = ikle      
   
  # {STRING} title
  # {OBJECT (name:str,unit:str)} var
  # {2D Array}
  # {2D Array(NELEM,3}
  def addMesh(self,title,var,points,ikle):
      self.empty = False
      self.addTITLE(title)
      self.addVAR(var)
      self.addPOIN(points)
      self.addIKLE(ikle)
              
  # {String}
  def writeSLF(self,output):
    if self.empty : self.createGrid()
    self.fole.update({ 'name': output })
    self.fole.update({ 'hook': open(output,'wb') })
    self.appendHeaderSLF()
    # ~~> Time stepping
    self.tags['times']=np.arange(self.values.shape[0])
    for t in range(self.NFRAME):
        self.appendCoreTimeSLF(t)
        self.appendCoreVarsSLF(self.values[t])
    self.fole['hook'].close()

  # {String}
  # {Float x6}
  def createGrid(self,title="Grid",xstart=-1,xend=1,xstep=0.1,ystart=-1,yend=1,ystep=0.1):
    xPoints = np.arange(xstart,xend+xstep,xstep)
    yPoints = np.arange(ystart,yend+ystep,ystep)
    
    xlen = len(xPoints)
    ylen = len(yPoints)
    
    xy = np.array([[x,y] for y in yPoints for x in xPoints],np.float32)
    
    ikle = []
    for row in range(ylen-1):
      for col in range(xlen-1):
        n1 = col+row*(ylen)
        n2 = (col+1)+row*(ylen)
        n3 = col+(row+1)*(ylen)
        n4 = (col+1)+(row+1)*(ylen)
        ikle.append([n1,n3,n2])
        ikle.append([n2,n3,n4])
    ikle =  np.array(ikle)
    self.addTITLE(title)
    self.addVAR({'name':'BOTTOM','unit':'m'})
    self.addPOIN(xy)
    self.addIKLE(ikle)
    self.tags['times']=np.arange(1)
    
  def printAtt(self):
    attr = {
      'NFRAME':self.NFRAME,
      'NVAR':self.NVAR,
      'NPOIN3':self.NPOIN3,
      'NELEM3':self.NELEM3,
      'EXTENT':self.EXTENT,
    }
    print(attr)