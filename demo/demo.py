import numpy as np
from slf import SLF

# Example 1
slf=SLF()                 # Create new Selafin file with a square grid
slf.writeSLF("demo1.slf") # Write new Selafin file

# Example 2
filename = 'demo1.slf'
slf=SLF(filename)         # Read Selafin file
slf.printAtt()            # Print attributes


# Example 3 - Change values
filename = 'demo1.slf'
slf=SLF(filename)
slf.values = slf.values + 1.0   # Change all values to 1.0
slf.writeSLF("demo2.0.slf")       # Write new Selafin file

for t in range(slf.NFRAME):
    for v in range(slf.NVAR):
        for p in range(slf.NPOIN3):
            slf.values[t,v,p] = (t*1000)+(v*500)+p
slf.writeSLF("demo2.1.slf")       # Write new Selafin file

for t in range(slf.NFRAME):
    for v in range(slf.NVAR):
        slf.values[t,v] =(t*1000)+(v*500) + (np.arange(slf.NPOIN3) / float(slf.NPOIN3-1))
slf.writeSLF("demo2.2.slf")       # Write new Selafin file

# Example 4 - Add variable
filename = 'demo1.slf'
slf=SLF(filename)
original = np.zeros((slf.NFRAME,slf.NVAR,slf.NPOIN3),np.float32)
for t in range(slf.NFRAME):
  original[t] = slf.getVALUES(t)

slf.addVAR({'name':'VARIABLE2','unit':'m'})
_values = np.arange(slf.NPOIN3) / float(slf.NPOIN3-1) # Add values
_values = _values[...,np.newaxis] + [0]             # Add var
_values = _values[...,np.newaxis]                     # Add time
_values = np.einsum('kji->ijk',_values)

slf.values = np.concatenate((original,_values), axis=1)
slf.writeSLF("demo3.slf")                             # Write new Selafin file
slf.printAtt()                                        # Print attributes

# Example 5 - Add time
filename = 'demo1.slf'
slf=SLF(filename)
original = np.zeros((slf.NFRAME,slf.NVAR,slf.NPOIN3),np.float32)
for t in range(slf.NFRAME):
  original[t] = slf.getVALUES(t)

values = np.zeros((10,slf.NVAR,slf.NPOIN3),np.float32)
for t in range(10):
    for ivar in range(slf.NVAR):
      values[t,ivar] = (np.arange(slf.NPOIN3) / float(slf.NPOIN3-1) + t + ivar)
slf.values = values


# Print
filename = 'demo1.slf'
slf=SLF(filename)         # Read Selafin file
slf.printInfo("demo1.js") # Print