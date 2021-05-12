# PySPuDS

PySPuDS is a python wrapper to [Structure Prediction and Diagnostic Software (SPuDS)](https://www.unf.edu/~michael.lufaso/spuds/) enabling high-throughput perovskite (Glazer tilt) structure predictions.

Currently PySPuDS only supports ABX<sub>3</sub> structure generation (SPuDS menu item 1). 
Future support for SPuDS menu items 3, 6, and 7 is planned.

## Installation
### SPuDS
Download and extract the [DOS, command line] version of [SPuDS](https://www.unf.edu/~michael.lufaso/spuds/) to your specified SPuDS-install-directory.
### PySPuDS
Download this repository (PySPuDS) and edit the SPuDS_install_directory variable in PySPuDS.py to point to your SPuDS-install-directory.

## Requirements
  - pymatgen

## Usage 
```python
from PySPuDS import SPuDS

A = {'Ca':2}
B = {'Ti':4}
X = {'O':-2}

Model = SPuDS(A,B,X,store_dir = 'PySPuDS_results')

for tilt in Model.allowed_tilts:   
    Model.write_default_input(tilt)
    Model.run()
    Model.write_cif()
    Model.store_results()
```
