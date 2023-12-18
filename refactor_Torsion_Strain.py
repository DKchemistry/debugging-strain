import os
import argparse
import csv
import refactor_TL_Functions as rtf

# Create the parser
parser = argparse.ArgumentParser(description="Debug SDF/Mol2 IO discrepency.")

# Add the arguments
parser.add_argument(
    "-i",
    "--InputFile",
    metavar="input",
    type=str,
    required=True,
    help="the input file (.mol2 or .sdf))",
)
parser.add_argument(
    "-o",
    "--OutputFile",
    metavar="output",
    type=str,
    required=True,
    help="the output file name (with extension)",
)

args = parser.parse_args()

# Check if the output file already exists
if os.path.exists(args.OutputFile):
    print("The output file already exists. Please choose a different name.")
    exit()

if args.InputFile[-5:] == ".mol2":
    names, ms = rtf.Mol2MolSupplier(args.InputFile)
elif args.InputFile[-4:] == ".sdf":
    names, ms = rtf.custom_sdmolsupplier_H(args.InputFile)
else:
    print("Error. Please pass a .mol2 or .sdf file.")
    exit()

print(str(len(names)) + " molecules finished reading. Calculating strain energy...")

with open(args.OutputFile, mode="w") as file_out:
    file_writer = csv.writer(file_out)
    count = 0
    for name in names:
        mol = ms[name]
        if mol is not None:
            try:
                M = rtf.TL_lookup(mol)
                mol_est = []
                mol_est += M.sum(0)
                mol_info = M.get_TPs()
                bond_info = [
                    item
                    for sublist in sorted(mol_info, key=lambda l: l[1], reverse=True)
                    for item in sublist
                ]
                file_writer.writerow([name] + mol_est + bond_info)
            except:
                count += 1
                file_writer.writerow([name] + ["NA"])
    file_out.close()

print(
    str(len(names) - count)
    + " successful / "
    + str(count)
    + " NA. Please check: "
    + args.OutputFile
)
