# Title
Description

## History
Empty

## Intent
Empty

## Installation
Empty

## Usage
```python
from slf import SLF
slf=SLF()                 # Create new Selafin file with a square grid
slf.writeSLF("demo1.slf") # Write new Selafin file
```

## Options
Empty

## Examples
### Print attributes
```python
from slf import SLF
slf=SLF()                 # Create new Selafin file with a square grid
slf.writeSLF("demo1.slf") # Write new Selafin file
slf.printAtt()
```
### Change values
```python
from slf import SLF
slf=SLF()                 # Create new Selafin file with a square grid
slf.values = slf.values + 1.0 
slf.writeSLF("demo1.slf") # Write new Selafin file
```
### Change values

## Test Cases
Empty

## License
ios-py is an open source python library and licensed under [MIT](../master/LICENSE).