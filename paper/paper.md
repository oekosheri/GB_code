---
title: 'GB\_code: A grain boundary generation code'
tags:
  - Python
  - grain boundary
  - crystallography
  - CSL
  - atomistic simulations (VASP, LAMMPS)
authors:
  - name: R. Hadian
    orcid: 0000-0002-9616-4602
    affiliation: 1
  - name: B. Grabowski
    orcid: 0000-0003-4281-5665
    affiliation: 1
  - name: J. Neugebauer
    affiliation: 1
affiliations:
 - name: Max-Planck-Institut fuer Eisenforschung, Duesseldorf, Germany
   index: 1
date: 19 July 2018
bibliography: paper.bib
---

# Summary

Grain boundaries (GBs) are crystalline borders between single crystals in materials microstructure. They play an important role in mechanical, chemical or electronic response of materials and are therefore essential to materials science and physics.

GBs are geometrical entities with a large parameter space that has been well formulated within a coincident site lattice (CSL) mathematical framework [@Sutton:1996]. One important computational advantage of the CSL formalism is that it enables the construction of GBs in a periodic setup for atomistic simulations. ``GB_code`` [@GB_code] uses the CSL construction to generate GB atomic structures (currently for cubic materials) systematically. It provides input atomic structures for large-scale atomistic simulations with interatomic potentials (as implemented e.g. in ``LAMMPS`` [@LAMMPS]) or _ab initio_, density-functional-theory (DFT) simulations (as implemented e.g. in ``VASP`` [@VASP]). These atomistic codes can further calculate different properties of the GBs. In addition to providing the input structures, the ``csl_generator.py`` script and the attached Jupyter notebooks have extra functionality to show how the CSL properties can be used to locate, classify and categorize different GBs and to extract detailed information about them.

``GB_code`` is mainly designed to be a command line tool as it is documented in detail in the README file of the repository
but the modules can also be accessed separately for example via the attached Jupyter notebooks. The code consists of two main scripts, ``csl_generator.py`` and ``gb_generator.py``, that should be used in this order to produce the final GB structures. The attached Jupyter notebooks in the Test directory, ``Usage_of_GB_code.ipynb`` and ``Dichromatic_pattern_CSL.ipynb``, input the two scripts as modules, the former addresses the general usage of the code with some extra tips and functions to locate GBs of interest, the latter depicts how CSL properties such as the overlapping patterns and displacement shift complete (DSC) vectors can be extracted and visualized. In the notebooks, two examples of the usage of the ``GB_code`` in our previous publications [@Pub1, @Pub2] have been shown as well.

``GB_code`` uses the analytical and mathematical formulations of the following works [@Sutton:1996, @Bollmann:1982, @Grimmer]. Some functionality from this code [@Marcin] on CSL has been used in a modified form in the ``GB_code``. 

### Statement of need:
``GB_code`` is an interactive toolbox to learn about grain boundaries and it is versatile for running high-throughput calculations. The target audience are students/scientists of materials science and physics at any level of familiarity with the topic. To our knowledge, due to the extensive usage of Python NumPy library, ``GB_code`` is faster than the existing python codes for GB generations. The user will be guided through in a step by step manner how to find and create the GB of interest. The code has been designed to be simple to use and instructive with a special attention to GB plane orientation which is often lacking in other grain boundary creation codes.

# Acknowledgements

R. Hadian would like to thank professor Mike Finnis for fruitful discussions. Funding from the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 639211) is also gratefully acknowledged.

# References
