
from dataclasses import dataclass, field
from typing import Literal
from rdkit.Chem import rdmolfiles, rdmolops, rdchem
import json
from icecream import ic
import subprocess


@dataclass
class ORCACalculator:
    name: str = "ORCA Calculator"
    base_name: str = "ORCA"
    executable: str = "orca"
    # Theory
    theory: str = "DFT"
    # Symmetry handling
    symmetry: bool = False
    # Initial Guess Options
    initial_guess: Literal["PATOM", "PMODEL", "HUECKEL", "HCORE"] = field(
        default="PATOM"
    )
    # Grid Options
    grid: Literal["DEFGRID1", "DEFGRID2", "DEFGRID3"] = field(default="DEFGRID2")
    # SCF Convergence
    scf_conv: Literal[
        "NORMALSCF",
        "LOOSESCF",
        "SLOPPYSCF",
        "STRONGSCF",
        "TIGHTSCF",
        "VERYTIGHTSCF",
        "EXTREMESCF",
    ] = field(
        default="NORMALSCF",
    )
    # Geometry Convergence
    geom_conv: Literal[
        "LOOSEOPT",
        "NORMALOPT",
        "TIGHTOPT",
        "VERYTIGHTOPT",
    ] = field(
        default="NORMALOPT",
    )
    # Convergence Acceleration
    conv_acc: Literal["DIIS", "KDIIS", "TRAH", "SOSCF", "DAMP", "LSHIFT"] = field(
        default="SOSCF",
    )
    CPCM: Literal[
        "1,1,1-trichloroethane",
        "1,1,2-trichloroethane",
        "1,2,4-trimethylbenzene",
        "1,2-dibromoethane",
        "1,2-dichloroethane",
        "1,2-ethanediol",
        "1,4-dioxane",
        "1-bromo-2-methylpropane",
        "1-bromooctane",
        "1-bromopentane",
        "1-bromopropane",
        "1-butanol",
        "1-chlorohexane",
        "1-chloropentane",
        "1-chloropropane",
        "1-decanol",
        "1-fluorooctane",
        "1-heptanol",
        "1-hexanol",
        "1-hexene",
        "1-hexyne",
        "1-iodobutane",
        "1-iodohexadecane",
        "1-iodopentane",
        "1-iodopropane",
        "1-nitropropane",
        "1-nonanol",
        "1-octanol",
        "1-pentanol",
        "1-pentene",
        "1-propanol",
        "2,2,2-trifluoroethanol",
        "2,2,4-trimethylpentane",
        "2,4-dimethylpentane",
        "2,4-dimethylpyridine",
        "2,6-dimethylpyridine",
        "2-bromopropane",
        "2-butanol",
        "2-chlorobutane",
        "2-heptanone",
        "2-hexanone",
        "2-methoxyethanol",
        "2-methyl-1-propanol",
        "2-methyl-2-propanol",
        "2-methylpentane",
        "2-methylpyridine",
        "2-nitropropane",
        "2-octanone",
        "2-pentanone",
        "2-propanol",
        "2-propen-1-ol",
        "e-2-pentene",
        "3-methylpyridine",
        "3-pentanone",
        "4-heptanone",
        "4-methyl-2-pentanone",
        "4-methylpyridine",
        "5-nonanone",
        "acetic acid",
        "acetone",
        "acetonitrile",
        "acetophenone",
        "ammonia",
        "aniline",
        "anisole",
        "benzaldehyde",
        "benzene",
        "benzonitrile",
        "benzyl alcohol",
        "bromobenzene",
        "bromoethane",
        "bromoform",
        "butanal",
        "butanoic acid",
        "butanone",
        "butanonitrile",
        "butyl ethanoate",
        "butylamine",
        "n-butylbenzene",
        "sec-butylbenzene",
        "tert-butylbenzene",
        "carbon disulfide",
        "carbon tetrachloride",
        "chlorobenzene",
        "chloroform",
        "a-chlorotoluene",
        "o-chlorotoluene",
        "conductor",
        "m-cresol",
        "o-cresol",
        "cyclohexane",
        "cyclohexanone",
        "cyclopentane",
        "cyclopentanol",
        "cyclopentanone",
        "decalin",
        "cis-decalin",
        "n-decane",
        "dibromomethane",
        "dibutylether",
        "o-dichlorobenzene",
        "e-1,2-dichloroethene",
        "z-1,2-dichloroethene",
        "dichloromethane",
        "diethyl ether",
        "diethyl sulfide",
        "diethylamine",
        "diiodomethane",
        "diisopropyl ether",
        "cis-1,2-dimethylcyclohexane",
        "dimethyl disulfide",
        "n,n-dimethylacetamide",
        "n,n-dimethylformamide",
        "dimethylsulfoxide",
        "diphenylether",
        "dipropylamine",
        "n-dodecane",
        "ethanethiol",
        "ethanol",
        "ethyl acetate",
        "ethyl methanoate",
        "ethyl phenyl ether",
        "ethylbenzene",
        "fluorobenzene",
        "formamide",
        "formic acid",
        "furan",
        "n-heptane",
        "n-hexadecane",
        "n-hexane",
        "hexanoic acid",
        "iodobenzene",
        "iodoethane",
        "iodomethane",
        "isopropylbenzene",
        "p-isopropyltoluene",
        "mesitylene",
        "methanol",
        "methyl benzoate",
        "methyl butanoate",
        "methyl ethanoate",
        "methyl methanoate",
        "methyl propanoate",
        "n-methylaniline",
        "methylcyclohexane",
        "n-methylformamide",
        "nitrobenzene",
        "nitroethane",
        "nitromethane",
        "o-nitrotoluene",
        "n-nonane",
        "n-octane",
        "n-pentadecane",
        "octanol(wet)",
        "pentanal",
        "n-pentane",
        "pentanoic acid",
        "pentyl ethanoate",
        "pentylamine",
        "perfluorobenzene",
        "phenol",
        "propanal",
        "propanoic acid",
        "propanonitrile",
        "propyl ethanoate",
        "propylamine",
        "pyridine",
        "tetrachloroethene",
        "tetrahydrofuran",
        "tetrahydrothiophene-s,s-dioxide",
        "tetralin",
        "thiophene",
        "thiophenol",
        "toluene",
        "trans-decalin",
        "tributylphosphate",
        "trichloroethene",
        "triethylamine",
        "n-undecane",
        "water",
        "xylene",
        "m-xylene",
        "o-xylene",
        "p-xylene",
        ""
    ] = field(
        default=None,
    )
    keepdens: bool = False
    keepints: bool = False
    density_functional: str = "B3LYP"
    basis_set: str = "6-31G"
    multiplicity: int = None
    charge: int = None
    parallel: Literal[
        "PAL2",
        "PAL3",
        "PAL4",
        "PAL5",
        "PAL6",
        "PAL7",
        "PAL8",
        "PAL16",
        "PAL32",
        "PAL64",
        ""
    ] = field(
        default=None,
    )
    aim: bool = False

    def get_settings(self):
        settings = {
            "theory": self.theory,
            "initial_guess": self.initial_guess,
            "runtype": self.runtype,
            "grid": self.grid,
            "scf_conv": self.scf_conv,
            "geom_conv": self.geom_conv,
            "conv_acc": self.conv_acc,
            "symmetry": self.symmetry,
            "keepdens": self.keepdens,
            "keepints": self.keepints,
        }
        return settings

    def get_properties(self, molecule: rdchem.Mol):
        with open(f"{self.base_name}.property.json", "r") as f:
            properties = json.load(f)
        geometry_keys = list(
            filter(lambda x: "Geometry" in x, properties.keys()),
        )
        final_geom = max(map(lambda x: int(x.split("_")[1]), geometry_keys))

        final_props = properties[f"Geometry_{final_geom}"]
        mayer_bond_orders = list(
            map(lambda x: x[0], final_props["Mayer_Population_Analysis"]["BONDORDERS"])
        )
        bonds = list(
            map(
                lambda x: [x[0], x[2]],
                final_props["Mayer_Population_Analysis"]["COMPONENTS"],
            )
        )
        bond_list = []
        for mbo, bond in zip(mayer_bond_orders, bonds):
            bond_list.append({"atom1": bond[0], "atom2": bond[1], "mbo": mbo})
        properties = {
            "global": {
                "Total Energy [eV]": final_props["SCF_Energy"]["SCF_ENERGY"]
                * 27.2113834,
                "Dipole Moment [D]": list(
                    map(
                        lambda x: x[0] * 2.541746473,
                        final_props["Dipole_Moment"]["DIPOLETOTAL"],
                    )
                ),
                "Dipole Magnitude [D]": final_props["Dipole_Moment"]["DIPOLEMAGNITUDE"]
                * 2.541746473,
            },
            "atomic": {
                "Mulliken Partial Charges [e]": list(
                    map(
                        lambda x: x[0],
                        final_props["Mulliken_Population_Analysis"]["ATOMICCHARGES"],
                    )
                ),
                "Loewdin Partial Charges [e]": list(
                    map(
                        lambda x: x[0],
                        final_props["Loewdin_Population_Analysis"]["ATOMICCHARGES"],
                    )
                ),
                "Mayer Partial Charges [e]": final_props["Mayer_Population_Analysis"][
                    "QA"
                ],
            },
            "bonds": {"Mayer Bond Order": bond_list},
        }
        return properties

    def run_orca(self, molecule: rdchem.Mol):
        keywords = [
            self.theory,
            self.initial_guess,
            self.runtype,
            self.grid,
            self.scf_conv,
            self.conv_acc,
        ]
        if self.runtype == "OPT":
            keywords.append(self.geom_conv)
        if self.symmetry:
            keywords.append("UseSym")
        if self.keepdens:
            keywords.append("KeepDens")
        if self.keepints:
            keywords.append("KeepInts")
        if self.parallel is not None:
            keywords.append(self.parallel)
        rdmolfiles.MolToXYZFile(molecule, "input.xyz")
        if self.charge is None:
            self.charge = rdmolops.GetFormalCharge(molecule)
        if self.multiplicity is None:
            total_electrons = 0
            if self.charge is not None:
                total_electrons += abs(self.charge)
            for atom in molecule.GetAtoms():
                total_electrons += atom.GetNumRadicalElectrons()
            total_electron_spin = total_electrons / 2
            multiplicity = int(2 * total_electron_spin + 1)
            if multiplicity == 0:
                multiplicity = 2
            self.multiplicity = multiplicity
        if self.aim:
            keywords.append("AIM")
        with open("input.inp", "w") as f:
            f.write(f"""%base "{self.base_name}"\n""")
            f.write(f"! {" ".join(keywords)}\n")
            f.write(f"""%Method\n\tWriteJSONPropertyfile True\nEnd\n""")
            f.write(f"* xyzfile {self.charge} {self.multiplicity} input.xyz\n")

        subprocess.call(f"{self.executable} input.inp > log.out", shell=True)
        return None
