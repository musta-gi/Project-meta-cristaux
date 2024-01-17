from Excel_Script_lattice import *


import pandas as pd
import shutil
from tempfile import NamedTemporaryFile

file_path = 'Generation_lattices.xlsx'

try:
    # Create a temporary file
    temp_file = NamedTemporaryFile(delete=False, suffix='.xlsx')
    temp_file_path = temp_file.name

    # Copy the content of the open file to the temporary file
    shutil.copyfile(file_path, temp_file_path)

    # Read the copied content using pandas
    df = pd.read_excel(temp_file_path)
    print(df)

    # Close and delete the temporary file
    temp_file.close()
except PermissionError as e:
    print("The file is open in another program and locked for editing:", e)
except FileNotFoundError:
    print("The file was not found.")
except Exception as e:
    print("An error occurred:", e)




# G7R
coordinates_G7R = 'G7R\coordonnees_motif_G7R.csv'
connectivity_G7R = 'G7R\connectivite_motif_G7R.csv'
# CFC
coordinates_CFC = 'CFC\coordonnees_motif_CFC.csv'
connectivity_CFC = 'CFC\connectivite_motif_CFC.csv'
# Octettruss
coordinates_Octettruss = 'Octettruss\coordonnees_motif_Octettruss.csv'
connectivity_Octettruss = 'Octettruss\connectivite_motif_Octettruss.csv'

coordinates, connectivity, name = 'Generation_lattices.xlsx', 'Generation_lattices.xlsx', 'CFC'#choose_structure('CFC')
#coordinates = coordinates_Octettruss
#connectivity = connectivity_Octettruss
constant = 6
width = 0.5
n = 3
nx, ny, nz = n, n, n  # Number of unit cells in each direction
#name = 'Octettruss'
destination = ''

axis = [1, 1, 1]
angle = np.pi/6

cells, beams = extract_data(coordinates, connectivity)
atoms, connectivity = generate_unit_cell(cells, beams, constant)


lattice, lattice_connectivity = duplicate_unit_cell(atoms, connectivity, constant, nx, ny, nz)

atoms, connectivity = lattice, lattice_connectivity #rotate_unit_cell(lattice, lattice_connectivity, axis, angle) #rotate the generated lattice
coordinates_DF, connectivity_DF, radius_DF = save_to_DF(atoms, connectivity, width)
number_of_struts = len(connectivity)
vol = number_of_struts*np.pi*(width**2)*constant
vol_glob = ((constant*n)**3)
density = vol/vol_glob
print(name)
print("number_of_struts", number_of_struts)
print("vol", vol)
print("vol_glob", vol_glob)
print("density", density)
#plot_lattice(atoms, connectivity) #<---- this is the one to plot
#plot_lattice(lattice, lattice_connectivity)
export_lattice(name, destination, coordinates_DF, connectivity_DF, radius_DF) #<---- this is the one to export
    