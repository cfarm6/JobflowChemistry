# JobflowChemistry

## Quantum Chemistry and Small Molecule Analysis within the Jobflow Framework

### Package Structure

To allow for introspection of the available tasks available in the package, the following submodules and jobs are available:

- Structure Inputs
  - PubChem CID Import
  - SMILES/SMARTS String
- Molecular Structure Generation
  - ETKDGv3 Generation
- Structure Energy Calculation
  - ML
    - AimNet2
    - MACE
    - OrbNet
  - DFT
    - ORCA
    - PySCF
  - xTB
    - GFN2-xTB
    - GFN1-xTB
    - GFN0-xTB
- Structure Properties
- Partial Charges
  - Multicharge

### Document Stores

The following set of document stores are implemented for managing data exchanges between jobs:

- files
  - System Structures:
    - .mol V3K files for small molecules
    - .sdf V3K for conformers included
    - .cif for X organic frameworks
    - .mmpdb for proteins
- settings
  - Job settings for mapping keywords args to values for reproducibility
- properties
  - global
    - Any properties relating to the overall system
  - atomic
    - Properties specific to the atoms involved
  - bond
    - Properties related to the concept of atomic bonding
