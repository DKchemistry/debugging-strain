from sys import argv #For passing commandline arguments
import csv #For writing a .csv file
# We will import the other required modules when we run TL_Functions.py

script, input = argv #Input from commandline argument

exec(open("TL_Functions.py").read()) #Run that script
if input[-5:] == ".mol2": #If fed in a .mol2 file
    names, ms = Mol2MolSupplier(input)
    # A list of names and a dictionary of mol2 objects with the names as keys
    # Make sure every mol2 object in the file has a different name, or not all
    # of them will end up in the final dictionary!
    output_name = input[:len(input)-5] #Everything but the ".mol2"
elif input[-4:] == ".db2": #If fed in a .db2 file
    exec(open("db2_to_mol2.py").read()) #Run that script
    file = db2_file_like(input)
    # Creates a string buffer object called "file" that looks like a .mol2 file
    names, ms = db2MolSupplier(file)
    output_name = input[:len(input)-4] #Everything but the ".db2"
elif input[-7:] == ".db2.gz": #If fed in a .db2.gz file
    exec(open("db2_to_mol2.py").read()) #Run that script
    file = db2_file_like(input)
    # Creates a string buffer object called "file" that looks like a .mol2 file
    names, ms = db2MolSupplier(file)
    output_name = input[:len(input)-7] #Everything but the ".db2.gz"
else:
    print("Error. Please pass a .mol2 or .db2 file.")
    exit()

print(str(len(names)) + " molecules finished reading. Calculating strain energy...")

# These (names and ms) will be our list of names and dictionary of molecules.
# We will loop over every molecule in this dictionary, create a TP_list object
# for each one, and output the information we want as a line in a .csv file.
# We will repeat the sum procedure for every cutoff value that we want
file_name = output_name + "_Torsion_Strain.csv"
with open(file_name, mode = "w") as file_out:
    # Initialize the output file
    file_writer = csv.writer(file_out) #Initialize writing with csv, using
    # https://realpython.com/python-csv/ as a guide
    count = 0
    for name in names: #Loop over every molecule in the list
        mol = ms[name] #Get the molecule with that name
        if mol is not None: #Check to make sure the molecule exists
            try:
                M = TL_lookup(mol) #Create a TP_list function
                mol_est = []
                mol_est += M.sum(0)
                mol_info = M.get_TPs() #The molecule's information
                bond_info = [item for sublist in sorted(mol_info, key = lambda l:l[1], reverse=True) for item in sublist]
                # Unlist the list of lists into just a list of elements
                file_writer.writerow([name] + mol_est + bond_info)
                # Add the list as a row
            except:
                count += 1
                file_writer.writerow([name] + ["NA"])
    file_out.close()

print(str(len(names)-count) + " successful / " + str(count) + " NA. Please check: " + file_name)

# END
