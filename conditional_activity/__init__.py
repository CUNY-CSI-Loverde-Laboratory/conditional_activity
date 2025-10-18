"""
Conditional_Activity
This code calculates the Conditional activity (time-correlated transitions) of a degree of freedom. In this case, the first sidechain dihedral angles (chi1)  for selected amino residues in protein (except ALA and GLY) and the base-sugar dihedral angle in DNA were used. The code seeks to find the kinetic correlation of amino acids side chains in the 3-dimensional native state of a protein and the base-sugar in a DNA strand. This module was built on MDAnalysis as a foundation using some functions in the MDAnalysis package.
"""

# Add imports here
from importlib.metadata import version

__version__ = version("conditional_activity")
