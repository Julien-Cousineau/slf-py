# SLFpy
Read and write selafin files in python

## Installation

This package was developed,tested and built using conda.
SLFpy uses numpy and matplotlib.
Only tested with python >=3.6

```bash
conda create -n slfpy python=3.8
conda activate slfpy
conda install -c meracan slfpy
```

Local installation
```bash
conda create -n slfpy python=3.8
conda activate slfpy

conda install -c conda-forge numpy 
git clone https://github.com/Julien-Cousineau/slf-py.git
pip install -e ./slf-py
```

## Usage
```python
from slfpy import SLF
# Create new Selafin file with a square grid
path = 'demo1.slf'
slf=SLF()                 
slf.write(path) 

# Read Selafin file and print attributes
slf=SLF(filename)         
slf.printAtt() 
```

### Usage, user guide and examples
[Docs](doc/doc_slfpy.ipynb)

### Testing
[Docs](test/README.md)

### License
[License](LICENSE)