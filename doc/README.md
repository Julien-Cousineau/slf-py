## Options
Empty

## Examples
### Change values
```python
slf.values = slf.values + 1.0

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
```
### Add values 
```python
values = np.zeros((10,slf.NVAR,slf.NPOIN3),np.float32)
for t in range(1,10):
    for ivar in range(slf.NVAR):
        values[t,ivar] = (np.arange(slf.NPOIN3) / float(slf.NPOIN3-1) + t + ivar)
slf.values = values

_values = np.arange(slf.NPOIN3) / float(slf.NPOIN3-1) # Add values
_values = values[...,np.newaxis] + [0,1]              # Add var
_values = values[...,np.newaxis]                      # Add time
_values = np.einsum('kij->ijk',values)


values = np.concatenate((original,rangec), axis=1)

```