# CONDACT
![CONDACT Logo](https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity/blob/main/CONDACT_Logo.png)

**CONDitional ACTivity: A Kinetic Framework for Mapping Allosteric Communication**

**Author:** Augustine C. Onyema\
**Assisted by:** Chukwuebuka Dikeocha\
**Lab:** Loverde Laboratory, CUNY\
**Status:** Validated & Accepted for Published 

--------------------------------------------------------------------------------------------------

## 1. Overview

CONDACT is an open-source Python library for quantifying **time-resolved
kinetic correlations** in biomolecular systems using molecular dynamics
(MD) simulations.

Unlike structural correlation methods (e.g., mutual information, PCA),
CONDACT focuses on **transition timing** of discrete geometric states
(e.g., side-chain dihedral angles) to uncover:

-   Residues with high **dynamical memory**
-   Directional **inter-residue kinetic coupling**
-   Long-range allosteric communication networks
-   Dynamically connected domains via principal eigenvectors
-   Statistically validated communication edges using permutation-based
    FDR

The method builds upon the concept of conditional activity and extends
it to large, heterogeneous systems including protein--DNA complexes.

The module is built as an MDKit. All files accepted by MDAnanlysis to create Universe can be used in CONDACT.\ 
Atom selection for MDAnalysis is also valid.

--------------------------------------------------------------------------------------------------

## 2. Core Concept: Conditional Activity

For two residues **X** and **Y**:

- **Transition times:** `T(X,i)`, `T(Y,i)`
- **Persistence time (memory):** `<tau_X>`
- **Exchange time:** `<tau_X|Y>`

**Conditional Activity:**

```
A[X][Y] = -log( <tau_X|Y> / <tau_X> )
```


For two residues X and Y:

-   Transition times: T(X,i), T(Y,i)\
-   Persistence time (memory): ⟨τ_X⟩\
-   Exchange time: ⟨τ_X\|Y⟩

Conditional Activity:

A\[X\]\[Y\] = -log( ⟨τ_X\|Y⟩ / ⟨τ_X⟩ )

Interpretation:

  Value              Meaning
  ------------------ ------------------------------------------
  A\[X\]\[Y\] \> 0   Transition in Y promotes transition in X
  A\[X\]\[Y\] = 0    Independent transitions
  A\[X\]\[Y\] \< 0   Suppression
  A\[X\]\[X\]        Dynamical memory of residue X

Importantly: A\[X\]\[Y\] ≠ A\[Y\]\[X\] (directional)

--------------------------------------------------------------------------------------------------

## 3. Repository Structure

```
conditional_activity/
├── analysis/
│   ├── condact.py        # Core CONDACT implementation
│   ├── fdr.py            # Permutation + BH-FDR correction
│   └── trace_route.py    # Communication pathway tracing
│
├── condact_tutorial/
│   └── tutorial.ipynb    # End-to-end example workflow
│
├── tests/                # Validation tests
├── pyproject.toml        # Build configuration
├── README.md             # This file
└── LICENSE
```

--------------------------------------------------------------------------------------------------

## 4. Installation

``` bash
git clone https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity.git
cd conditional_activity
pip install .
```

Dependencies: - Python ≥ 3.9 - MDAnalysis - NumPy - SciPy - Pandas - Matplotlib - Jupyter Notebook

--------------------------------------------------------------------------------------------------

## 5. Basic Workflow

### Step 1 --- Load Trajectory

``` python
import MDAnalysis as mda
from conditional_activity.analysis.condact import CONDACT

u = mda.Universe("topology.prmtop", "trajectory.xtc")
```

### Step 2 --- Run Conditional Activity

``` python
ca = CONDACT(u,
            selected_resid='1-5',               # Select the residues of interest eg. here 1 to 3
            No_of_peaks_protein=3,              # The number of peaks for amino acid as integer
            peak_boundaries_protein=[120, 240], # list of boundaries for amino acid
            states_protein="XYZ",               # States for amino acids in str
            No_of_peaks_nucleic=3,              # The number of peaks for nucleotides if present as integer
            peak_boundaries_nucleic=[140, 330], # Slist of boundaries for nucleotides if present
            states_nucleic="SAS",               # States for nucleic acid
            saving_frequency=10,                # Saving frequency of trajectory in picoseconds in int
            keep_negative=False                 # All negative dihedral will be converted to secondary angles
            )
CA_matrix = ca.mean_conditional_activity()
```
Read the work & Visit tutorial.ipynb to learn how to assign states and analyze trajectory

### Step 3 --- Statistical Validation (FDR)

``` python
from conditional_activity.analysis.fdr import CONDACT_FDR

fdr = CONDACT_FDR(CA_matrix, n_perm=1000, q=0.05)
results = fdr.run()
```

--------------------------------------------------------------------------------------------------

## 6. Advanced Analyses

### Dynamical Memory Mapping

Diagonal entries: A\[X\]\[X\]

### Inter-residue conditional activity Mapping

Offdiagonal enteries: A\[X\]\[Y\]

### Principal Eigenvector Analysis

Predicts dynamically connected domains.

--------------------------------------------------------------------------------------------------

## 7. Why CONDACT?

✔ Directional\
✔ Kinetically grounded\
✔ Long-range sensitive\
✔ Statistically validated\
✔ Works for protein & DNA\
✔ Predicts communication hubs\
✔ Mutation-sensitive

--------------------------------------------------------------------------------------------------

## 8. Limitations

✔ Reliable conditional activity estimates require residues to sample multiple state transitions; sparse transitions can reduce statistical stability. Our module adopted at least 10 transitions.
✔ Discretizing continuous dihedral angles into finite states introduces coarse-graining that may overlook subtle conformational variations.
✔ The quadratic scaling of residue–residue comparisons increases computational cost for very large biomolecular systems.

--------------------------------------------------------------------------------------------------

## 9. Citation

If you use CONDACT, please cite:

Onyema et al., *Mapping Allosteric Communication in the Nucleosome with
Conditional Activity*.

------------------------------------------------------------------------

## 10. Data Availability

-   GitHub:
    https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity
-   Trajectories: https://zenodo.org/records/17193099
-   GNU General Public License, version 2 (see the file [LICENSE](https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity/blob/main/LICENSE)).


#### Copyright

Copyright (c) 2025, Augustine Onyema, Loverde Laboratory, CUNY

### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.
Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using Conditional_Activity in published work.