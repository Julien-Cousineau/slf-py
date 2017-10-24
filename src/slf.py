#! /usr/bin/python
import os,sys
import numpy as np
from struct import unpack,pack

class SLF(object):
  def __init__(self,filePath=''):
    self.empty = True
    self.filePath = filePath
    self._extent = None
    self._trixy = None
    self._triarea = None
    self._tricentroid = None
    self._values = None     
    
    self.initialized(filePath)

  def initialized(self,fileName):
    self.file = {}
    self.file.update({ 'name': fileName })
    self.file.update({ 'endian': ">" })    # "<" means little-endian, ">" means big-endian
    self.file.update({ 'float': ('f',4) }) #'f' size 4, 'd' = size 8
    if fileName != '':
       self.file.update({ 'hook': open(fileName,'rb') })
       # ~~> checks endian encoding
       self.file['endian'] = self.getEndianFromChar(self.file['hook'],80)
       # ~~> header parameters
       self.tags = { 'meta': self.file['hook'].tell() } #TODO remove?
       self.getHeaderMetaDataSLF()
       # ~~> sizes and connectivity
       self.getHeaderIntegersSLF()
       # ~~> checks float encoding
       self.file['float'] = self.getFloatTypeFromFloat(self.file['hook'],self.file['endian'],self.NPOIN3)
       # ~~> xy mesh
       self.getHeaderFloatsSLF()
       # ~~> time series
       self.tags = { 'cores':[],'times':[] }
       self.getTimeHistorySLF()
    else:
       self.TITLE = ''
       self.NBV1 = 0; self.NBV2 = 0; self.NVAR = self.NBV1 + self.NBV2
       self.VARINDEX = range(self.NVAR)
       self.IPARAM = []
       self.NELEM3 = 0; self.NPOIN3 = 0; self.NDP3 = 0; self.NPLAN = 1
       self.NELEM2 = 0; self.NPOIN2 = 0; self.NDP2 = 0
       self.NBV1 = 0; self.VARNAMES = []; self.VARUNITS = []
       self.NBV2 = 0; self.CLDNAMES = []; self.CLDUNITS = []
       self.IKLE3 = []; self.IKLE2 = []; self.IPOB2 = []; self.IPOB3 = []; self.MESHX = []; self.MESHY = []
       self.tags = { 'cores':[],'times':[] }
    self.fole = {}
    self.fole.update({ 'name': '' })
    self.fole.update({ 'endian': self.file['endian'] })
    self.fole.update({ 'float': self.file['float'] })
    self.tree = None
    self.neighbours = None
    self.edges = None
    self.alterZnames = []

  def getEndianFromChar(self,f,nchar):
    pointer = f.tell()
    endian = ">"       # "<" means little-endian, ">" means big-endian
    l,c,chk = unpack(endian+'i'+str(nchar)+'si',f.read(4+nchar+4))
    if chk!=nchar:
      endian = "<"
      f.seek(pointer)
      l,c,chk = unpack(endian+'i'+str(nchar)+'si',f.read(4+nchar+4))
    if l!=chk:
      print '... Cannot read '+str(nchar)+' characters from your binary file'
      print '     +> Maybe it is the wrong file format ?'
      sys.exit(1)
    f.seek(pointer)
    return endian

  def getFloatTypeFromFloat(self,f,endian,nfloat):
    pointer = f.tell()
    ifloat = 4
    cfloat = 'f'
    l = unpack(endian+'i',f.read(4))
    if l[0]!=ifloat*nfloat:
      ifloat = 8
      cfloat = 'd'
    r = unpack(endian+str(nfloat)+cfloat,f.read(ifloat*nfloat))
    chk = unpack(endian+'i',f.read(4))
    if l!=chk:
      print '... Cannot read '+str(nfloat)+' floats from your binary file'
      print '     +> Maybe it is the wrong file format ?'
      sys.exit(1)
    f.seek(pointer)
    return cfloat,ifloat

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Parsing the Big- and Little-Endian binary file
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def getHeaderMetaDataSLF(self):
    f = self.file['hook']
    endian = self.file['endian']
    # ~~ Read title ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l,self.TITLE,chk = unpack(endian+'i80si',f.read(4+80+4))
    # ~~ Read NBV(1) and NBV(2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l,self.NBV1,self.NBV2,chk = unpack(endian+'iiii',f.read(4+8+4))
    self.NVAR = self.NBV1 + self.NBV2
    self.VARINDEX = range(self.NVAR)
    # ~~ Read variable names and units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    self.VARNAMES = []; self.VARUNITS = []
    for _ in range(self.NBV1):
       l,vn,vu,chk = unpack(endian+'i16s16si',f.read(4+16+16+4))
       self.VARNAMES.append(vn)
       self.VARUNITS.append(vu)
    self.CLDNAMES = []; self.CLDUNITS = []
    for _ in range(self.NBV2):
       l,vn,vu,chk = unpack(endian+'i16s16si',f.read(4+16+16+4))
       self.CLDNAMES.append(vn)
       self.CLDUNITS.append(vu)
    # ~~ Read IPARAM array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    d = unpack(endian+'12i',f.read(4+40+4))
    self.IPARAM = np.asarray( d[1:11] )
    # ~~ Read DATE/TIME array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    self.DATETIME = [1972,07,13,17,15,13]
    if self.IPARAM[9] == 1:
       d = unpack(endian+'8i',f.read(4+24+4))
       self.DATETIME = np.asarray( d[1:9] )

  def getHeaderIntegersSLF(self):
    f = self.file['hook']
    endian = self.file['endian']
    # ~~ Read NELEM3, NPOIN3, NDP3, NPLAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l,self.NELEM3,self.NPOIN3,self.NDP3,self.NPLAN,chk = unpack(endian+'6i',f.read(4+16+4))
    self.NELEM2 = self.NELEM3
    self.NPOIN2 = self.NPOIN3
    self.NDP2 = self.NDP3
    self.NPLAN = max( 1,self.NPLAN )
    if self.IPARAM[6] > 1:
       self.NPLAN = self.IPARAM[6] # /!\ How strange is that ?
       self.NELEM2 = self.NELEM3 / ( self.NPLAN - 1 )
       self.NPOIN2 = self.NPOIN3 / self.NPLAN
       self.NDP2 = self.NDP3 / 2
    # ~~ Read the IKLE array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.seek(4,1)
    self.IKLE3 = np.array( unpack(endian+str(self.NELEM3*self.NDP3)+'i',f.read(4*self.NELEM3*self.NDP3)) ) - 1
    f.seek(4,1)
    self.IKLE3 = self.IKLE3.reshape((self.NELEM3,self.NDP3))
    if self.NPLAN > 1: self.IKLE2 = np.compress( np.repeat([True,False],self.NDP2), self.IKLE3[0:self.NELEM2], axis=1 )
    else: self.IKLE2 = self.IKLE3
    # ~~ Read the IPOBO array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.seek(4,1)
    self.IPOB3 = np.asarray( unpack(endian+str(self.NPOIN3)+'i',f.read(4*self.NPOIN3)) )
    f.seek(4,1)
    self.IPOB2 = self.IPOB3[0:self.NPOIN2]

  def getHeaderFloatsSLF(self):
    f = self.file['hook']
    endian = self.file['endian']
    # ~~ Read the x-coordinates of the nodes ~~~~~~~~~~~~~~~~~~
    ftype,fsize = self.file['float']
    f.seek(4,1)
    self.MESHX = np.asarray( unpack(endian+str(self.NPOIN3)+ftype,f.read(fsize*self.NPOIN3))[0:self.NPOIN2] )
    f.seek(4,1)
    # ~~ Read the y-coordinates of the nodes ~~~~~~~~~~~~~~~~~~
    f.seek(4,1)
    self.MESHY = np.asarray( unpack(endian+str(self.NPOIN3)+ftype,f.read(fsize*self.NPOIN3))[0:self.NPOIN2] )
    f.seek(4,1)

  def getTimeHistorySLF(self):
    f = self.file['hook']
    endian = self.file['endian']
    ftype,fsize = self.file['float']
    ATs = []; ATt = []
    while True:
       try:
          ATt.append(f.tell())
          # ~~ Read AT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          f.seek(4,1)
          ATs.append(unpack(endian+ftype,f.read(fsize))[0])
          f.seek(4,1)
          # ~~ Skip Values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          f.seek(self.NVAR*(4+fsize*self.NPOIN3+4),1)
       except:
          ATt.pop(len(ATt)-1)   # since the last record failed the try
          break
    self.tags.update({ 'cores': ATt })
    self.tags.update({ 'times': np.asarray(ATs) })

  def getVariablesAt( self,frame,varsIndexes ):
    f = self.file['hook']
    endian = self.file['endian']
    ftype,fsize = self.file['float']
    if fsize == 4: z = np.zeros((len(varsIndexes),self.NPOIN3),dtype=np.float32)
    else: z = np.zeros((len(varsIndexes),self.NPOIN3),dtype=np.float64)
    # if tags has 31 frames, len(tags)=31 from 0 to 30, then frame should be >= 0 and < len(tags)
    if frame < len(self.tags['cores']) and frame >= 0:
       f.seek(self.tags['cores'][frame])
       f.seek(4+fsize+4,1)
       for ivar in range(self.NVAR):
          f.seek(4,1)
          if ivar in varsIndexes:
             z[varsIndexes.index(ivar)] = unpack(endian+str(self.NPOIN3)+ftype,f.read(fsize*self.NPOIN3))
          else:
             f.seek(fsize*self.NPOIN3,1)
          f.seek(4,1)
    return z

  def alterEndian(self):
    if self.fole['endian'] == ">": self.fole['endian'] = "<"
    else: self.fole['endian'] = ">"
  
  def alterFloat(self):
    if self.fole['float'] == ('f',4): self.fole['float'] = ('d',8)
    else: self.fole['float'] = ('f',4)
  
  def alterVALUES(self,vars=None,mZ=1,pZ=0):
    if vars != None:
       self.alterZm = mZ; self.alterZp = pZ; self.alterZnames = vars.split(':')

  def appendHeaderSLF(self):
    f = self.fole['hook']
    endian = self.fole['endian']
    ftype,fsize = self.fole['float']
    # ~~ Write title ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'i80si',80,self.TITLE,80))
   # ~~ Write NBV(1) and NBV(2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'iiii',4+4,self.NBV1,self.NBV2,4+4))
    # ~~ Write variable names and units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i in range(self.NBV1):
       f.write(pack(endian+'i',32))
       f.write(pack(endian+'16s',self.VARNAMES[i]))
       f.write(pack(endian+'16s',self.VARUNITS[i]))
       f.write(pack(endian+'i',32))
    for i in range(self.NBV2):
       f.write(pack(endian+'i',32))
       f.write(pack(endian+'16s',self.CLDNAMES[i]))
       f.write(pack(endian+'16s',self.CLDUNITS[i]))
       f.write(pack(endian+'i',32))
    # ~~ Write IPARAM array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'i',4*10))
    for i in range(len(self.IPARAM)): f.write(pack(endian+'i',self.IPARAM[i]))
    f.write(pack(endian+'i',4*10))
    # ~~ Write DATE/TIME array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if self.IPARAM[9] == 1:
       f.write(pack(endian+'i',4*6))
       for i in range(6): f.write(pack(endian+'i',self.DATETIME[i]))
       f.write(pack(endian+'i',4*6))
    # ~~ Write NELEM3, NPOIN3, NDP3, NPLAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'6i',4*4,self.NELEM3,self.NPOIN3,self.NDP3,1,4*4))  #/!\ where is NPLAN ?
    # ~~ Write the IKLE array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'i',4*self.NELEM3*self.NDP3))
    f.write(pack(endian+str(self.NELEM3*self.NDP3)+'i',*(self.IKLE3.ravel()+1)))
    f.write(pack(endian+'i',4*self.NELEM3*self.NDP3))
    # ~~ Write the IPOBO array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'i',4*self.NPOIN3))
    f.write(pack(endian+str(self.NPOIN3)+'i',*(self.IPOB3)))
    f.write(pack(endian+'i',4*self.NPOIN3))
    # ~~ Write the x-coordinates of the nodes ~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'i',fsize*self.NPOIN3))
    #f.write(pack(endian+str(self.NPOIN3)+ftype,*(np.tile(self.MESHX,self.NPLAN))))
    for i in range(self.NPLAN): f.write(pack(endian+str(self.NPOIN2)+ftype,*(self.MESHX)))
    f.write(pack(endian+'i',fsize*self.NPOIN3))
    # ~~ Write the y-coordinates of the nodes ~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'i',fsize*self.NPOIN3))
    #f.write(pack(endian+str(self.NPOIN3)+ftype,*(np.tile(self.MESHY,self.NPLAN))))
    for i in range(self.NPLAN): f.write(pack(endian+str(self.NPOIN2)+ftype,*(self.MESHY)))
    f.write(pack(endian+'i',fsize*self.NPOIN3))
  
  def appendCoreTimeSLF( self,t ):
    f = self.fole['hook']
    endian = self.fole['endian']
    ftype,fsize = self.fole['float']
    # Print time record
    if type(t) == type(0.0): f.write(pack(endian+'i'+ftype+'i',fsize,t,fsize))
    else: f.write(pack(endian+'i'+ftype+'i',fsize,self.tags['times'][t],fsize))

  def appendCoreVarsSLF( self,VARSOR ):
    f = self.fole['hook']
    endian = self.fole['endian']
    ftype,fsize = self.fole['float']
    # Print variable records
    for v in VARSOR:
       f.write(pack(endian+'i',fsize*self.NPOIN3))
       f.write(pack(endian+str(self.NPOIN3)+ftype,*(v)))
       f.write(pack(endian+'i',fsize*self.NPOIN3))
  
  def putContent( self,fileName):
    self.fole.update({ 'name': fileName })
    self.fole.update({ 'hook': open(fileName,'wb') })
    ibar = 0
    self.appendHeaderSLF()
    for t in range(len(self.tags['times'])):
       ibar += 1
       self.appendCoreTimeSLF(t)
       self.appendCoreVarsSLF(self.getVALUES(t))
    self.fole['hook'].close()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Tool Box
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def getVALUES( self,t ):
    VARSOR = self.getVariablesAt( t,self.VARINDEX )
    for v in self.alterZnames:
       for iv in range(len(self.VARNAMES)):
          if v.lower() in self.VARNAMES[iv].lower(): VARSOR[iv] = self.alterZm * VARSOR[iv] + self.alterZp
       for iv in range(len(self.CLDNAMES)):
          if v.lower() in self.CLDNAMES[iv].lower(): VARSOR[iv+self.NBV1] = self.alterZm * VARSOR[iv+self.NBV1] + self.alterZp
    return VARSOR
  
  def getSERIES( self,nodes,varsIndexes=[]):
    f = self.file['hook']
    endian = self.file['endian']
    ftype,fsize = self.file['float']
    if varsIndexes == []: varsIndexes = self.VARINDEX
    # ~~ Ordering the nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This assumes that nodes starts at 1
    onodes = np.sort(np.array( zip(range(len(nodes)),nodes), dtype=[ ('0',int),('1',int) ] ),order='1')
    # ~~ Extract time profiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if fsize == 4: z = np.zeros((len(varsIndexes),len(nodes),len(self.tags['cores'])),dtype=np.float32)
    else: z = np.zeros((len(varsIndexes),len(nodes),len(self.tags['cores'])),dtype=np.float64)
    f.seek(self.tags['cores'][0])
    for t in range(len(self.tags['cores'])):
       f.seek(self.tags['cores'][t])
       f.seek(4+fsize+4,1)
       for ivar in range(self.NVAR):
          f.seek(4,1)
          if ivar in varsIndexes:
             jnod = onodes[0]
             f.seek(fsize*(jnod[1]-1),1)
             z[varsIndexes.index(ivar),jnod[0],t] = unpack(endian+ftype,f.read(fsize))[0]
             for inod in onodes[1:]:
                f.seek(fsize*(inod[1]-jnod[1]-1),1)
                z[varsIndexes.index(ivar),inod[0],t] = unpack(endian+ftype,f.read(fsize))[0]
                jnod = inod
             f.seek(fsize*self.NPOIN3-fsize*jnod[1],1)
          else:
             f.seek(fsize*self.NPOIN3,1)
          f.seek(4,1)
    return z

  def __del__(self):
    if self.file['name'] != '': self.file['hook'].close()

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