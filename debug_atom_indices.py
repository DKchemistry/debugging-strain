 
from refactor_TL_Functions import Mol2MolSupplier, custom_sdmolsupplier_H

def compare_molecules(mol1, mol2):
    atoms1 = [(atom.GetIdx(), atom.GetSymbol()) for atom in mol1.GetAtoms()]
    atoms2 = [(atom.GetIdx(), atom.GetSymbol()) for atom in mol2.GetAtoms()]

    if atoms1 != atoms2:
        print(f"Discrepancy found: {atoms1} vs {atoms2}")
        return False

    return True

# Load the molecules with your custom function
names_custom, mols_custom = custom_sdmolsupplier_H("output.sdf")

# Load the molecules with Mol2MolSupplier
names_mol2, mols_mol2 = Mol2MolSupplier("example.mol2")

# Iterate over all molecules
for name in names_custom:
    mol_custom = mols_custom[name]
    mol_mol2 = mols_mol2[name]

    print(f"Molecule: {name}")
    if not compare_molecules(mol_custom, mol_mol2):
        print("Molecules do not match!")
    else:
        print("Molecules match!")