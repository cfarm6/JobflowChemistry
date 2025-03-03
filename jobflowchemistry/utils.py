from rdkit.Chem import rdchem, rdDetermineBonds
from ase import Atoms
import re
import json


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
    for i, pos in enumerate(positions):
        conf.SetAtomPosition(i, pos)
    return mol


def parse_gaussian_output(file_path):
    # Initialize data structures
    frequencies = []
    reduced_masses = []
    force_constants = []
    ir_intensities = []
    raman_activity = []
    depolarization = []
    displacement_vectors = {}

    # Open and read the Gaussian output file
    with open(file_path, "r") as f:
        lines = f.readlines()

    # Flags for identifying relevant sections
    in_reduced_masses = False
    in_force_constants = False
    in_ir_intensities = False
    in_raman_activity = False
    in_depolarization = False
    in_displacement_vectors = False
    current_mode = 0

    # Process each line of the file
    for line in lines:
        # Check for frequency section
        if "Frequencies --" in line:
            in_force_constants = False
            in_ir_intensities = False
            in_raman_activity = False
            in_depolarization = False
            in_displacement_vectors = False

            # continue
            freqs = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line)))
            frequencies.extend(freqs)
            in_reduced_masses = True
            continue

        # Extract reduced masses
        if in_reduced_masses:
            if "Red. masses --" in line:
                reduced_masses.extend(
                    list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line)))
                )
                in_reduced_masses = False
                in_force_constants = True
            continue
        # Extract Force constants
        if in_force_constants:
            if "Frc consts" in line:
                force_constants.extend(
                    list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line)))
                )
                in_force_constants = False
                in_ir_intensities = True
            continue

        # Extract IR intensities
        if in_ir_intensities:
            if "IR Inten" in line:
                ir_intensities.extend(
                    list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line)))
                )
                in_ir_intensities = False
                in_raman_activity = True
            continue

        # Extract Raman Activity
        if in_raman_activity:
            if "Raman Activ" in line:
                raman_activity.extend(
                    list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line)))
                )
                in_raman_activity = False
                in_depolarization = True
            continue

        # Extract Depolar
        if in_depolarization:
            if "Depolar" in line:
                depolarization.extend(
                    list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line)))
                )
                in_depolarization = False
                in_displacement_vectors = True
            continue
        if line.lstrip().startswith("Atom"):
            in_displacement_vectors = True
            # Skip the header line
            continue
        if in_displacement_vectors:
            if re.match(r"^\s*\d", line):  # Look for lines with atom numbers
                displacement_data = list(
                    map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line))
                )
                if len(displacement_data) < 5:
                    # When the next section is reached, reset flags
                    in_displacement_vectors = False
                    current_mode += len(x)
                    continue
                atom_id = int(displacement_data[0])  # The first entry is the atom ID
                element_type = int(displacement_data[1])
                xyz = displacement_data[2:]
                x = xyz[0::3]
                y = xyz[1::3]
                z = xyz[2::3]
                for i in range(len(x)):
                    if str(current_mode + i) not in displacement_vectors:
                        displacement_vectors[str(current_mode + i)] = {}
                    for xyz in zip(xyz):
                        displacement_vectors[str(current_mode + i)][str(atom_id)] = {
                            "x": x[i],
                            "y": y[i],
                            "z": z[i],
                        }
    return (
        frequencies,
        reduced_masses,
        force_constants,
        ir_intensities,
        raman_activity,
        depolarization,
        displacement_vectors,
    )


def fix_json_file(input_file, output_file):
    try:
        with open(input_file, "r") as file:
            json_data = file.read()

        # Remove trailing commas in arrays and objects
        json_data = re.sub(r",\s*([\]}])", r"\1", json_data)

        # Remove leading commas in arrays or objects (if any malformed structures)
        json_data = re.sub(r"([{\[])\s*,", r"\1", json_data)

        # Fix specific issues like an empty array with a trailing comma
        json_data = json_data.replace(",\n\n   ", ",\n   ")

        # Fix specific malformed array entries (e.g., the "bond orders" issue)
        json_data = re.sub(r"\[\s*,", "[", json_data)

        # Attempt to parse the cleaned string as JSON
        data = json.loads(json_data)

        # Write the cleaned and parsed JSON data back to a file with proper formatting
        with open(output_file, "w") as file:
            json.dump(data, file, indent=3)

        print(f"JSON file fixed and saved to: {output_file}")

    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
