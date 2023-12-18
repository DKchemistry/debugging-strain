import refactor_TL_Functions_debug as rtf

# Load the molecules with your custom function
names_custom, mols_custom = rtf.custom_sdmolsupplier_H("output.sdf")

# Load the molecules with Mol2MolSupplier
names_mol2, mols_mol2 = rtf.Mol2MolSupplier("example.mol2")


def compare_lists(list1, list2, name, property_name):
    # Check if the lengths of the outer lists are the same
    if len(list1) != len(list2):
        print(f"Different number of elements in {property_name} for molecule {name}")
    else:
        # Iterate over the lists and compare each element
        difference_found = False
        for i in range(len(list1)):
            if list1[i] != list2[i]:
                print(
                    f"Difference found in element {i} of {property_name} for molecule {name}"
                )
                print(f"TL_customSDF: {list1[i]}, TL_mol2: {list2[i]}")
                difference_found = True
        if not difference_found:
            print(f"No differences found in {property_name} for molecule {name}")


for name in names_custom:
    mol_custom = mols_custom[name]
    mol_mol2 = mols_mol2[name]

    print(f"Molecule in sdf: {name}")
    TL_custom = rtf.TL_lookup(mol_custom)
    print(f"Molecule in mol2: {name}")
    TL_mol2 = rtf.TL_lookup(mol_mol2)

    compare_lists(TL_custom.indeces, TL_mol2.indeces, name, "indeces")
    #compare_lists(TL_custom.angles, TL_mol2.angles, name, "angles")
    compare_lists(TL_custom.smarts, TL_mol2.smarts, name, "smarts")
    #compare_lists(TL_custom.hc, TL_mol2.hc, name, "hc")
    #compare_lists(TL_custom.methods, TL_mol2.methods, name, "methods")
    #compare_lists(TL_custom.E, TL_mol2.E, name, "E")
    #compare_lists(TL_custom.CI_l, TL_mol2.CI_l, name, "CI_l")
    #compare_lists(TL_custom.CI_u, TL_mol2.CI_u, name, "CI_u")
    #compare_lists(TL_custom.flags, TL_mol2.flags, name, "flags")
