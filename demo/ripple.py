from slf import SLF
import numpy as np

class mySLF(SLF):
    def __init__(self,filePath=''):
        SLF.__init__(self, filePath)
        
    def AddRipple(self,frames=np.array([0])):
        x = self.MESHX[...,np.newaxis]
        y = self.MESHY[...,np.newaxis]
        
        frames = frames / float(len(frames)) * 4.* np.pi # 2 cycles
        values = np.sin(2.*(np.power(x,2.)+np.power(y,2.)+frames))/2.
        values = values[...,np.newaxis]
        values = np.einsum('ijk->jki', values)
        self.values = values


if __name__ == "__main__":
    slf=mySLF()
    slf.AddRipple()
    slf.writeSLF("../data/ripple.2D.1.slf")
    slf.AddRipple(np.arange(100))
    slf.writeSLF("../data/ripple.2D.100.slf")