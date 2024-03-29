# PySPuDS

**PySPuDS** is a python wrapper to [Structure Prediction and Diagnostic Software](https://www.unf.edu/~michael.lufaso/spuds/) (SPuDS) enabling high-throughput perovskite (Glazer tilt) structure predictions.

Currently **PySPuDS** only supports ABX<sub>3</sub> structure generation, i.e. SPuDS menu item 1. 

run_spuds_v1.py has limited support for menu items 1, 3, 6, and 7. This version is static with the release of the [Theoretical Multinary Perovskite Oxides Dataset](https://contribs.materialsproject.org/projects/Multinary_Oxides)



# Installation
### SPuDS
Download and extract the [DOS, command line] version of [SPuDS](https://www.unf.edu/~michael.lufaso/spuds/) to your specified SPuDS-install-directory.

Note: SPuDS is only supported on Windows os

### PySPuDS
Download this repository (**PySPuDS**) and copy ```PySPuDS.py``` and ```ABX3_SPuDS_symops.json``` to your specified SPuDS-install-directory.



# Requirements
  - Windows os
  - SPuDS (DOS version > 2.21.05.11)
  - pymatgen




# Usage
SPuDS requires only A/B site assignments, elements, and oxidation states to predict Glazer tilt structures. Similarly, the ```SPuDS()``` class imported from ```PySPuDS.py``` requires only dictionaries of A, B, and X elements with corresponding oxidation states. The ```store_dir``` variable defaults to 'PySPuDS_results', but can be changed upon initialization of the ```SPuDS()``` class.


For a given Glazer tilt system, **PySPuDS**: 1) writes an input file (defaults to no Jahn-Teller distortions at 298 K ), 2) runs the SPuDS program using the generated input file, 3) creates a pymatgen-compatible .cif structure file using symmetry operations specified in ```ABX3_SPuDS_symops.json```, and 4) stores the results in the specified ```store_dir``` using a default naming scheme.

# run_spuds_v1.py
This version has similar usage just replace the import line with:
```python
from run_spuds_v1 import SPuDS
```

### Example Usage 
```python
from PySPuDS import SPuDS

A = {'Ca':2} #dict of cation A {element: oxidation state}
B = {'Ti':4} #dict of cation B {element: oxidation state}
X = {'O':-2} #dict of anion X {element: oxidation state}

store_directory = 'PySPuDS_results' #where to store SPuDS output and .cif files

Model = SPuDS(A, B, X, store_dir = store_directory)

for tilt in Model.allowed_tilts:   
    Model.write_default_input(tilt)
    Model.run()
    Model.write_cif()
    Model.store_results()
```
