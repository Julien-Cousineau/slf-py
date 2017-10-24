import numpy as np
from slf import SLF
slf=SLF()                 # Create new Selafin file with a square grid
slf.writeSLF("demo1.slf") # Write new Selafin file

filename = 'demo1.slf'
slf=SLF(filename)         # Create new Selafin file with a square grid
slf.writeSLF("demo2.slf") # Write new Selafin file
# slf.printAtt()            # Print attributes

slf.values = slf.values + 1.0   # Change all values to 1.0
slf.writeSLF("demo3.slf")       # Write new Selafin file

for t in range(slf.NFRAME):
    for v in range(slf.NVAR):
        for p in range(slf.NPOIN3):
            slf.values[t,v,p] = (t*1000)+(v*500)+p

for t in range(slf.NFRAME):
    for v in range(slf.NVAR):
        slf.values[t,v] =(t*1000)+(v*500) + (np.arange(slf.NPOIN3) / float(slf.NPOIN3-1))

slf.addVAR({'name':'VARIABLE2','unit':'m'})
_values = np.arange(slf.NPOIN3) / float(slf.NPOIN3-1) # Add values
_values = _values[...,np.newaxis] + [0,1]             # Add var
_values = _values[...,np.newaxis]                     # Add time
slf.values = np.einsum('kij->jik',_values)            # Change the shape (frame,var,points)
print slf.values
slf.writeSLF("demo3.slf")       # Write new Selafin file
slf.printAtt()            # Print attributes