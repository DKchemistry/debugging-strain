import refactor_TL_Functions_debug_2 as rtf2

mol2_names, mol2_mols = rtf2.Mol2MolSupplier("example.mol2")
# make variable one_mol that is equal to the first value of the first key in mol2_mols
one_mol = mol2_mols[list(mol2_mols.keys())[0]]
print(one_mol)
# Call TL_lookup on one_mol
rtf2.TL_lookup(one_mol)

