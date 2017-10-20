from slf import SLF
import numpy as np

def createEmptySLF(filename):
  slf=SLF()
  slf.writeSLF(filename)
  slf.printAtt()

def readSLF(filename):
  slf = SLF(filename)
  slf.printAtt()
  
def changeValuesV1(input,output):
  slf = SLF(input)
  slf.values = slf.values + 1.0 
  slf.writeSLF(output)
  slf.printAtt()

def changeFrame(input,output):
  slf = SLF(input)
  
  range = np.arange(slf.NPOIN3) / float(slf.NPOIN3-1)
  values = range[...,np.newaxis] + [0,1]
  values = values[...,np.newaxis]
  values = np.einsum('kij->ijk',values)
  slf.values =values
  slf.writeSLF(output)
  slf.printAtt()

def addFrames(input,output):
  slf = SLF(input)
  
  range = np.arange(slf.NPOIN3) / float(slf.NPOIN3-1)
  values = range[...,np.newaxis] + [0,1]
  values = values[...,np.newaxis]
  values = np.einsum('kij->ijk',values)
  slf.values = values
  slf.writeSLF(output)
  slf.printAtt()

def addVars(input,output):
  slf = SLF(input)
  
  slf.addVAR({'name':'VARIABLE2','unit':'m'})
  
  range = np.arange(slf.NPOIN3) / float(slf.NPOIN3-1)
  values = range[...,np.newaxis] + [0,1]
  values = values[...,np.newaxis]
  values = np.einsum('kij->jik',values)
  slf.values = values
  slf.writeSLF(output)
  slf.printAtt()


def addVarToOriginal(input,output):
  slf = SLF(input)
  original = np.zeros((slf.NFRAME,slf.NVAR,slf.NPOIN3),np.float32)
  for t in range(slf.NFRAME):
      original[t] = slf.getVALUES(t)
      
  slf.addVAR({'name':'area','unit':'m2'})
  
  rangea = np.arange(slf.NPOIN3) / float(slf.NPOIN3-1)
  
  ranget = rangea[np.newaxis,np.newaxis,...] + [[[0]],[[1]]]
  
  rangeb = rangea[...,np.newaxis] + [0,1]
  rangec = rangeb[...,np.newaxis]
  rangec = np.einsum('kij->ijk',rangec)
  

  values = np.concatenate((original,rangec), axis=1)

  slf.values = values
  slf.writeSLF(output)
  slf.printAtt()

def addFrameVars(input,output):
  slf = SLF(input)
  
  slf.addVAR({'name':'VARIABLE2','unit':'m'})
  values = np.zeros((10,slf.NVAR,slf.NPOIN3),np.float32)
  for t in range(10):
    for ivar in range(slf.NVAR):
      values[t,ivar] = (np.arange(slf.NPOIN3) / float(slf.NPOIN3-1) + t + ivar)
  
  slf.values = values
  slf.writeSLF(output)
  print("readSLF - NPOIN3:{0},NVAR:{1}, NELEM3:{2}, shape:{3}".format(slf.NPOIN3,slf.NVAR,slf.NELEM3,slf.values.shape))
  print("addVars - min:{0}, max:{1}".format(np.min(slf.values),np.max(slf.values)))

if __name__ == "__main__":
  filename = "demo1.slf"
  filename2 = "demo2.slf"
  filename3 = "demo3.slf"
  filename4 = "demo4.slf"
  filename5 = "demo5.slf"
  filename6 = "demo6.slf"
  createSLF(filename)
  # readSLF(filename)
  changeFrame(filename,filename2)
  # addFrames(filename2,filename3)
  # addVars(filename2, filename4)
 
  # addFrameVars(filename2, filename5)
  addVarToOriginal(filename2, filename6)
  # readSLF(filename2)

  