from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools

from qupkake import predict as qpka
from qupkake.mol_dataset import MolDataset
import os
# 0. Setup folders
if not os.path.exists("raw"):
    os.mkdir("raw")
if not os.path.exists("processed"):
    os.mkdir("processed")
smiles = "C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F"
# 1. smiles to sdf
# 1.a embed molecule
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)
mol.SetProp("_Name", "molecule")
with Chem.SDWriter("raw/mol.sdf") as f:
    f.write(mol)
# 1.b set name
# 1.c update args
filename = "mol.sdf"
mol_col = "ROMol"
name_col = "name"
root = "."
output = "output.sdf"
mp = False
tautomerize=False
# 2. Convert the sdf to a mol_dataset
dataset = MolDataset(
    root=root,
    filename=filename,
    tautomerize=tautomerize,
    name_col=name_col,
    mol_col=mol_col,
    mp=mp,
)
prot_model, deprot_model, pka_model = qpka.load_models()
prot_indices = qpka.predict_sites(dataset, prot_model)
deprot_indices = qpka.predict_sites(dataset, deprot_model)
qpka.make_sites_prediction_files(root, dataset, prot_indices, deprot_indices, output)
pair_dataset = qpka.load_mol_pair_dataset(
    root=root,
    filename=output,
    name_col=name_col,
    mol_col="ROMol",
    idx_col="idx",
    type_col="pka_type",
    mp=mp,
)
pka_predictions = qpka.predict_pka(pair_dataset, pka_model)
print(pka_predictions.numpy())
