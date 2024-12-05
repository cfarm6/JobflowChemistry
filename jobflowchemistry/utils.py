from rdkit.Chem import rdchem, rdDetermineBonds
from ase import Atoms

def rdkit2ase(molecule: rdchem.Mol) -> Atoms:
    atomic_numbers = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    positions = molecule.GetConformer().GetPositions()
    return Atoms(numbers=atomic_numbers, positions=positions)

def ase2rdkit(atoms: Atoms, mol: rdchem.Mol = None):
    if mol is None:
        atomic_numbers = atoms.get_atomic_numbers()
        rw_mol = rdchem.RWMol()
        for atomic_number in atomic_numbers:
            rw_mol.AddAtom(rdchem.Atom(int(atomic_number)))
        mol = rdchem.Mol(rw_mol)
        mol.AddConformer(rdchem.Conformer(len(atomic_numbers)), assignId=True)
    positions = atoms.get_positions()
    conf = mol.GetConformer()
    for i,pos in enumerate(positions):
        conf.SetAtomPosition(i, pos)
    return mol
