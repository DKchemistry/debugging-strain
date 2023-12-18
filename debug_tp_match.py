import refactor_TL_Functions_debug as rtf

# Load the molecules with your custom function
names_custom, mols_custom = rtf.custom_sdmolsupplier_H("output.sdf")

# Load the molecules with Mol2MolSupplier
names_mol2, mols_mol2 = rtf.Mol2MolSupplier("example.mol2")

for name in names_custom:
    mol_custom = mols_custom[name]
    mol_mol2 = mols_mol2[name]

    print(f"Molecule in sdf: {name}")
    TL_custom = rtf.TL_lookup(mol_custom)

    # Store the contents of match_debug_list in a variable
    match_debug_list_custom = list(rtf.match_debug_list)

    # Clear match_debug_list at the end of each loop iteration
    rtf.match_debug_list.clear()

    print(f"Molecule in mol2: {name}")
    TL_mol2 = rtf.TL_lookup(mol_mol2)

    # Compare the elements of match_debug_list_custom and match_debug_list
    for custom, mol2 in zip(match_debug_list_custom, rtf.match_debug_list):
        if custom != mol2:
            print(f"Difference found: {custom} != {mol2}")

    # Clear match_debug_list at the end of each loop iteration
    rtf.match_debug_list.clear()