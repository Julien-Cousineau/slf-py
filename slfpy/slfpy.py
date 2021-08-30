#! /usr/bin/python
import os,sys
import numpy as np
from struct import unpack,pack
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt

class SLF(object):
  def __init__(self,filePath=''):
    self.filePath = filePath
    self._extent = None
    self._trixy = None
    self._triarea = None
    self._tricentroid = None
    self._values = None
    self._kdtree = None
    
    self.initialized(filePath)

  def initialized(self,fileName):
    self.file = {}
    self.file.update({ 'name': fileName })
    self.file.update({ 'endian': ">" })    # "<" means little-endian, ">" means big-endian
    self.file.update({ 'float': ('f',4) }) #'f' size 4, 'd' = size 8
    if fileName != '':
       self.empty = False
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
       self.empty = True
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
      # self._createGrid()
    self.fole = {}
    self.fole.update({ 'name': '' })
    self.fole.update({ 'endian': self.file['endian'] })
    self.fole.update({ 'float': self.file['float'] })
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
      print('... Cannot read '+str(nchar)+' characters from your binary file')
      print('     +> Maybe it is the wrong file format ?')
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
      print('... Cannot read '+str(nfloat)+' floats from your binary file')
      print('     +> Maybe it is the wrong file format ?')
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
    self.TITLE = self.TITLE.decode('ascii').strip()
    # ~~ Read NBV(1) and NBV(2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    l,self.NBV1,self.NBV2,chk = unpack(endian+'iiii',f.read(4+8+4))
    self.NVAR = self.NBV1 + self.NBV2
    self.VARINDEX = range(self.NVAR)
    # ~~ Read variable names and units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    self.VARNAMES = []; self.VARUNITS = []
    for _ in range(self.NBV1):
       l,vn,vu,chk = unpack(endian+'i16s16si',f.read(4+16+16+4))
       vn = vn.decode('ascii').strip()
       vu = vu.decode('ascii').strip()
       self.VARNAMES.append(vn)
       self.VARUNITS.append(vu)
    self.CLDNAMES = []; self.CLDUNITS = []
    for _ in range(self.NBV2):
       l,vn,vu,chk = unpack(endian+'i16s16si',f.read(4+16+16+4))
       vn = str(vn)
       vu = str(vu)
       self.CLDNAMES.append(vn)
       self.CLDUNITS.append(vu)
    # ~~ Read IPARAM array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    d = unpack(endian+'12i',f.read(4+40+4))
    self.IPARAM = np.asarray( d[1:11] )
    # ~~ Read DATE/TIME array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    self.DATETIME = [1972,7,13,17,15,13]
    if self.IPARAM[9] == 1:
       d = unpack(endian+'8i',f.read(4+24+4))
       print(d)
       self.DATETIME = np.asarray( d[1:7] )
       

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
       self.NELEM2 = int(self.NELEM3 / ( self.NPLAN - 1 ))
       self.NPOIN2 = int(self.NPOIN3 / self.NPLAN)
       self.NDP2 = int(self.NDP3 / 2)
    # ~~ Read the IKLE array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.seek(4,1)
    self.IKLE3 = np.array( unpack(endian+str(self.NELEM3*self.NDP3)+'i',f.read(4*self.NELEM3*self.NDP3)) ) - 1
    f.seek(4,1)
    self.IKLE3 = self.IKLE3.reshape((self.NELEM3,self.NDP3))
    # print(self.NELEM2)
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
    self.MESHX = np.asarray( unpack(endian+str(self.NPOIN3)+ftype,f.read(fsize*self.NPOIN3))[0:self.NPOIN2],np.float64)
    f.seek(4,1)
    # ~~ Read the y-coordinates of the nodes ~~~~~~~~~~~~~~~~~~
    f.seek(4,1)
    self.MESHY = np.asarray( unpack(endian+str(self.NPOIN3)+ftype,f.read(fsize*self.NPOIN3))[0:self.NPOIN2],np.float64)
    f.seek(4,1)

  def dt2cal(self,dt):
    """
    Convert array of datetime64 to a calendar array of year, month, day, hour,
    minute, seconds, microsecond with these quantites indexed on the last axis.

    Parameters
    ----------
    dt : datetime64 array (...)
        numpy.ndarray of datetimes of arbitrary shape

    Returns
    -------
    cal : uint32 array (..., 7)
        calendar array with last axis representing year, month, day, hour,
        minute, second, microsecond
    """
  
    # allocate output
    out = np.empty(dt.shape + (7,), dtype="u4")
    # decompose calendar floors
    Y, M, D, h, m, s = [dt.astype("M8[{0}]".format(x)) for x in "YMDhms"]
    out[..., 0] = Y + 1970  # Gregorian Year
    out[..., 1] = (M - Y) + 1  # month
    out[..., 2] = (D - M) + 1  # dat
    out[..., 3] = (dt - D).astype("m8[h]")  # hour
    out[..., 4] = (dt - h).astype("m8[m]")  # minute
    out[..., 5] = (dt - m).astype("m8[s]")  # second
    out[..., 6] = (dt - s).astype("m8[us]")  # microsecond
    return out
  
  def setDatetime(self,value):
    self.DATETIME = self.dt2cal(value)[:6]
    self.IPARAM[9] = 1
  
  def getDatetime(self):
    print(self.DATETIME)
    [y,M,d,h,m,s]=self.DATETIME
    return np.datetime64('{0:04d}-{1:02d}-{2:02d}T{3:02d}:{4:02d}:{5:02d}'.format(y,M,d,h,m,s)) + self.tags['times'].astype('timedelta64[s]')
  
    

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
    z = np.zeros((len(varsIndexes),self.NPOIN3),self.dtype)
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

  def getVALUES( self,t ):
    VARSOR = self.getVariablesAt( t,self.VARINDEX )
    for v in self.alterZnames:
       for iv in range(len(self.VARNAMES)):
          if v.lower() in self.VARNAMES[iv].lower(): VARSOR[iv] = self.alterZm * VARSOR[iv] + self.alterZp
       for iv in range(len(self.CLDNAMES)):
          if v.lower() in self.CLDNAMES[iv].lower(): VARSOR[iv+self.NBV1] = self.alterZm * VARSOR[iv+self.NBV1] + self.alterZp
    return VARSOR

  @property
  def XY(self):
    return np.column_stack((self.MESHX,self.MESHY))

  @property
  def kdtree(self):
    import time
    from scipy import spatial
    if self._kdtree is None:
      self._kdtree = spatial.KDTree(self.XY)
      
    return self._kdtree
    
  def closestIndexes(self,pts):
    return self.kdtree.query(pts)[1]

  
  def getSERIES( self,nodes,varsIndexes=[]):
    f = self.file['hook']
    endian = self.file['endian']
    ftype,fsize = self.file['float']
    if varsIndexes == []: varsIndexes = self.VARINDEX
    # ~~ Ordering the nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This assumes that nodes starts at 1
    onodes=np.sort(np.array(list(zip(range(len(nodes)), nodes)), dtype=[('0', int), ('1', int)]), order='1')
    # onodes = np.sort(np.array( zip(range(len(nodes)),nodes), dtype=[ ('0',int),('1',int) ] ),order='1')
    # ~~ Extract time profiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    z = np.zeros((len(varsIndexes),len(nodes),len(self.tags['cores'])),self.dtype)
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
    # print(self.VARNAMES)
    f.write(pack(endian+'i80si',80,self.TITLE.encode(),80))
   # ~~ Write NBV(1) and NBV(2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'iiii',4+4,self.NBV1,self.NBV2,4+4))
    # ~~ Write variable names and units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i in range(self.NBV1):
       f.write(pack(endian+'i',32))
       f.write(pack(endian+'16s',self.VARNAMES[i].encode()))
       f.write(pack(endian+'16s',self.VARUNITS[i].encode()))
       f.write(pack(endian+'i',32))
    for i in range(self.NBV2):
       f.write(pack(endian+'i',32))
       f.write(pack(endian+'16s',self.CLDNAMES[i].encode()))
       f.write(pack(endian+'16s',self.CLDUNITS[i].encode()))
       f.write(pack(endian+'i',32))
    # ~~ Write IPARAM array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.write(pack(endian+'i',4*10))
    # print(self.IPARAM)
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
    # print(self.IKLE3.ravel())
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
    # print(self.tags['times'][t])

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


  def __del__(self):
    if self.file['name'] != '': self.file['hook'].close()

  def getVarsIndexes(self,names):
    return list(np.where(np.char.strip(np.asarray(self.VARNAMES)) == np.asarray(names)[:, np.newaxis])[1])

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
      values = np.zeros((self.NFRAME,self.NVAR,self.NPOIN3),self.dtype)
      if not self.empty:
        for t in range(self.NFRAME):
          values[t] = self.getVALUES(t)
      self.empty = False
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
    if len(self.tags['times'])==0:self.tags['times']=np.arange(1)
    return len(self.tags['times'])
  
  @property
  def dtype(self):
    ftype,fsize = self.file['float']
    dtype = np.float32
    if (fsize == 8):dtype = np.float64
    return dtype
  
  @property
  def TRIXY(self):
    if self._trixy is None:
        x,y=(self.MESHX,self.MESHY)
        z = np.zeros(x.size,np.float64)
        indexes = self.IKLE2
        self._trixy = np.concatenate((x[indexes][...,np.newaxis],y[indexes][...,np.newaxis],z[indexes][...,np.newaxis]), axis=2)
    return self._trixy
  
  @property
  def TRIAREA(self):
    if self._triarea is None:
      trixy = self.TRIXY
      v0 = trixy[:,1,] - trixy[:,0,:]
      v1 = trixy[:,2,:] - trixy[:,0,:]
      self._triarea = 0.5* np.cross(v0,v1)[:,2]
    return self._triarea

  @property
  def TRICENTROID(self):
    if self._tricentroid is None:
      self._tricentroid = np.mean(self.TRIXY,axis=1)
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
  
  def addVARS(self, vars):
    for var in vars:
      self.addVAR(var)
    
  # {2D Array}
  def addPOIN(self,points):
    self.IPOB2 = np.arange(len(points),dtype="i4")
    self.IPOB3 = self.IPOB2
    self.IPARAM = np.zeros(10,np.int32)
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
    self.IKLE2 = ikle.astype("int")
    self.IKLE3 = ikle.astype("int")      
   
  # {STRING} title
  # {OBJECT (name:str,unit:str)} var
  # {2D Array}
  # {2D Array(NELEM,3}
  def addMesh(self,points,ikle,**kwargs):
      title=kwargs.get("title","Mesh")
      var =kwargs.get("var",{"name":"BOTTOM","unit":"m"})
      values =kwargs.get("values",None)
      ipobo =kwargs.get("ipobo",None)
      self.empty = False
      self.addTITLE(title)
      if isinstance(var,list):
        for v in var: self.addVAR(v)
      else:self.addVAR(var)
      self.addPOIN(points)
      self.addIKLE(ikle)
      if values is not None:
        self.values=values
      else:
        self.values = np.zeros((self.NFRAME,self.NVAR,self.NPOIN3),self.dtype)
      
      if ipobo is not None:
        self.IPOB3=ipobo
      # else:
      #   self.IPOB3=np.zeros(self.NPOIN2)
      
      return self
      
    # {String}
  def writeSLFperFrame(self, output,func):
    self.fole.update({'name': output})
    self.fole.update({'hook': open(output, 'wb')})
    self.appendHeaderSLF()
    func(self)
    self.fole['hook'].close()

    # {String}
  def writeHeader(self, output):
      self.fole.update({'name': output})
      self.fole.update({'hook': open(output, 'wb')})
      self.appendHeaderSLF()
  
  def writeFrame(self,time,frame):
    self.appendCoreTimeSLF(time)
    self.appendCoreVarsSLF(frame)
      
  # {String}
  def write(self,output):
    self.fole.update({ 'name': output })
    self.fole.update({ 'hook': open(output,'wb') })
    self.appendHeaderSLF()
    # ~~> Time stepping
    self.tags['times']=np.arange(self.values.shape[0])
    for t in range(self.NFRAME):
        self.appendCoreTimeSLF(t)
        self.appendCoreVarsSLF(self.values[t])
    self.fole['hook'].close()
  
  
  @staticmethod
  def createGrid(xstart=-1,xend=1,xstep=0.1,ystart=-1,yend=1,ystep=0.1):
    xPoints = np.arange(xstart,xend+xstep,xstep)
    yPoints = np.arange(ystart,yend+ystep,ystep)
    
    xlen = len(xPoints)
    ylen = len(yPoints)
  
    x, y = np.meshgrid(xPoints, yPoints)
    x=x.ravel()
    y=y.ravel()
    elem = []
    for row in range(ylen-1):
      for col in range(xlen-1):
        n1 = col+row*(xlen)
        n2 = (col+1)+row*(xlen)
        n3 = col+(row+1)*(xlen)
        n4 = (col+1)+(row+1)*(xlen)
        elem.append([n1,n3,n2])
        elem.append([n2,n3,n4])
    elem =  np.array(elem)  
    return {"x":x,"y":y,"elem":elem}
  
  @staticmethod
  def plotGrid(x,y,elem,filePath):
    tri = Triangulation(x, y, elem.astype("int32"))
    plt.triplot(tri)
    plt.savefig(filePath)
    
  @staticmethod
  def createGridOld(title="Grid",xstart=-1,xend=1,xstep=0.1,ystart=-1,yend=1,ystep=0.1):
    xPoints = np.arange(xstart,xend+xstep,xstep)
    yPoints = np.arange(ystart,yend+ystep,ystep)
    
    xlen = len(xPoints)
    ylen = len(yPoints)
    
    xy = np.array([[x,y] for y in yPoints for x in xPoints],np.float64)
    
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
    return {"title":title,"xy":xy,"ikle":ikle}
  
  def _createGrid(self,*args,**kwargs):
    obj=SLF.createGridOld()
    self.addTITLE(obj['title'])
    self.addVAR({'name':'BOTTOM','unit':'m'})
    self.addPOIN(obj['xy'])
    self.addIKLE(obj['ikle'])
    self.tags['times']=np.arange(1)
    self.values # get values
    self.empty = False
    
  def printAtt(self):
    attr = {
      'NFRAME':self.NFRAME,
      'NVAR':self.NVAR,
      'NPOIN3':self.NPOIN3,
      'NELEM3':self.NELEM3,
      'EXTENT':self.EXTENT,
    }
    print(attr)
    print(self.VARNAMES)
    
    
    
  def printInfo(self,filename):
    ftype,fsize = self.file['float']
    # form = {'float':'{:12.6e}'.format} if (ftype=='f') else {'float':'{:20.14e}'.format}
    # np.set_printoptions(threshold=np.nan,formatter=form)
    # np.set_printoptions(threshold=np.nan,suppress=True)
    with open(filename, 'w') as f:
      f.write("exports.NELEM3=%d;\n" % self.NELEM3)
      f.write("exports.NPOIN3=%d;\n" % self.NPOIN3)
      f.write("exports.NFRAME=%d;\n" % self.NFRAME)
      f.write("exports.NFRAME10=%d;\n" % 10)
      f.write("exports.MESHX=new Float32Array({});\n".format(np.array2string(self.MESHX.flatten(),separator=', ')))
      f.write("exports.MESHY=new Float32Array({});\n".format(np.array2string(self.MESHY.flatten(),separator=', ')))
      f.write("exports.TRIXY=new Float32Array({});\n".format(np.array2string(self.TRIXY.flatten(),separator=', ')))
      f.write("exports.ELEMENTS=new Uint32Array({});\n".format(np.array2string(self.IKLE3.flatten(),separator=', ')))
      f.write("exports.TRIAREA=new Float32Array({});\n".format(np.array2string(self.TRIAREA.flatten(),separator=', ')))
      f.write("exports.CX=new Float32Array({});\n".format(np.array2string(self.TRICENTROID[:,0].flatten(),separator=', ')))
      f.write("exports.CY=new Float32Array({});\n".format(np.array2string(self.TRICENTROID[:,1].flatten(),separator=', ')))
      
      frame10 = (np.arange(self.NPOIN3) / float(self.NPOIN3)) + 9
      f.write("exports.FRAME10=new Float32Array({});\n".format(np.array2string(frame10.flatten(),separator=', ')))