from slfpy import SLF
import numpy as np
def test_slfpy():
    slf1=SLF()                 # Create new Selafin file with a square grid
    slf1.writeSLF("demo1.slf") # Write new Selafin file
    slf2=SLF("demo1.slf") # Write new Selafin file
    
    np.testing.assert_almost_equal(slf1.MESHX,slf2.MESHX,decimal=6)
    np.testing.assert_almost_equal(slf1.MESHY,slf2.MESHY,decimal=6)
    np.testing.assert_almost_equal(slf1.IKLE3,slf2.IKLE3)
    
    obj=SLF.createGrid()
    np.testing.assert_almost_equal(slf1.MESHX,obj['xy'][:,0],decimal=6)
    np.testing.assert_almost_equal(slf1.MESHY,obj['xy'][:,1],decimal=6)
    np.testing.assert_almost_equal(slf1.IKLE3,obj['ikle'])
    
    
if __name__ == "__main__":
  test_slfpy()