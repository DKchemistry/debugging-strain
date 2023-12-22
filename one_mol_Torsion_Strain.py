import refactor_TL_Functions_debug_2 as rtf2

mol2_names, mol2_mols = rtf2.Mol2MolSupplier("example.mol2")
# make variable one_mol that is equal to the first value of the first key in mol2_mols
one_mol = mol2_mols[list(mol2_mols.keys())[0]]
two_mol = mol2_mols[list(mol2_mols.keys())[1]]
three_mol = mol2_mols[list(mol2_mols.keys())[2]]
# Call TL_lookup on one_mol
rtf2.TL_lookup(one_mol)

print(f"Total XML SMARTS processed: {rtf2.global_XML_smarts_counter}")

# Reset the counter
rtf2.global_XML_smarts_counter = 0

# rtf2.TL_lookup(two_mol)

# print(f"Total XML SMARTS processed: {rtf2.global_XML_smarts_counter}")
