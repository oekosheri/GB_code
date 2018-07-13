# GB_code
This python package helps you create orthogonal grain boundary supercells for atomistic calculations. The code is based on the 
coincident site lattice (CSL) formulations for cubic materials (sc, bcc, fcc, diamond). I intend to extend it to hcp structures soon.

# Structure
There are two main scripts: [_csl_generator_](./csl_generator.py) and [_gb_generator_](./csl_generator.py) which you need to use in this order to produce the final gb structure.
In this description I will explain the steps to use the code in the Terminal and I have also attached two _jupyter notebooks_ ([Usage_of_GB_code.ipynb](./Usage_of_GB_code.ipynb), [Dichromatic_pattern_CSL_.ipynb](./Dichromatic_pattern_CSL_.ipynb)) which
describe how the code be can be accessed and used in the notebboks by various examples. 

# Usage
To pick a grain boundary (GB) 5 degrees of freedom need to be fixed: rotation axis, rotation angle and GB plane orientation.
We start by choosing only an axis, say a low index [1, 1, 1] and list the possibilities for the angle (sigma). Once you pick
this a CSL minimal cell will be created and you can see a list of possible GB plane orientations within the CSL that you pick
from. In the jupyter notebooks, [Usage_of_GB_code.ipynb](./Usage_of_GB_code.ipynb), example criteria have been shown to help pinpoint 
the boundary plane of interest.









