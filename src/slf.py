from fromOpenTelemac.parserSELAFIN import *

class SLF(SELAFIN):
    def __init__(self,filePath=''):
        self.filePath = filePath
        self.slfExist()
        
        SELAFIN.__init__(self, filePath)
        self._area = None
        self._centroid = None
        self._values = None        
    