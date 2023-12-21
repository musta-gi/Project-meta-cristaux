import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

connectivite_cube = pd.read_csv('Naked cube G7R\connectivite_motif_nu.csv', header=None)
#convert to list of tuples
connectivite_cube = [tuple(x) for x in connectivite_cube.values]
print(connectivite_cube)

corners_cube = pd.read_csv('coordonnees_cube.csv')

#convert to a ndarray
corners_cube = np.array(corners_cube)

def generate_fcc_unit_cell(a):
    """
    Generate a face-centered cubic (FCC) unit cell.

    Parameters:
    a (float): Lattice constant, the length of the unit cell cube.

    Returns:
    numpy.ndarray: Coordinates of atoms in the unit cell.
    """
    # Coordinates for the corners and face centers of the cube
    '''corners = a * np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                            [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])
    face_centers = a * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                                 [0.5, 0.5, 1], [0.5, 1, 0.5], [1, 0.5, 0.5]])
    '''
    # coordinates for the FCC
    '''corners = a * np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0],
                            [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1]])
    face_centers = a * np.array([[0.5, 0, 0.5], [0, 0.5, 0.5], [1, 0.5, 0.5], [0.5, 1, 0.5],
                                 [0.5, 0.5, 1], [0.5, 0.5, 0]])'''
    # coordinates for the G7R
    corners = a * np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0],
                            [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1]])
    face_centers = a * np.array([[0.5, 0.5, 0.5]])

    atoms = np.vstack((corners, face_centers))
    print("atoms", atoms) 

    # Define connectivity within the unit cell
    connectivity = []
    connectivity = connectivite_cube
    '''for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if np.linalg.norm(atoms[i] - atoms[j]) < 1.1 * a / np.sqrt(2):  # threshold for nearest neighbors
                connectivity.append((i, j))'''
    
    print("atoms", atoms)
    print("connectivity", connectivity)
    return atoms, connectivity

def duplicate_unit_cell(unit_cell, connectivity, nx, ny, nz):
    """
    Duplicate the unit cell to form a lattice with connectivity.

    Parameters:
    unit_cell (numpy.ndarray): Coordinates of atoms in the unit cell.
    connectivity (list): List of tuples representing connected atom indices.
    nx, ny, nz (int): Number of unit cells in the x, y, and z directions.

    Returns:
    numpy.ndarray: Coordinates of atoms in the lattice.
    list: Connectivity of the lattice.
    """
    lattice = []
    lattice_connectivity = []
    num_atoms_per_cell = len(unit_cell)

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                # Translation vector for the current unit cell
                translation = np.array([x, y, z]) * a
                # Translate and add the unit cell to the lattice
                translated_cell = unit_cell + translation
                lattice.append(translated_cell)

                # Add connectivity for this unit cell
                offset = num_atoms_per_cell * (x * ny * nz + y * nz + z)
                for conn in connectivity:
                    lattice_connectivity.append((conn[0] + offset, conn[1] + offset))
    #print("connectivity", lattice_connectivity)
    
    return np.vstack(lattice), lattice_connectivity
    

# Lattice constant (you can adjust this)
a = 6 # desnité = a/(2*r)

# Generate the unit cell and its connectivity
unit_cell, unit_cell_connectivity = generate_fcc_unit_cell(a)

# Duplicate the unit cell to create a lattice (change nx, ny, nz as needed)
nx, ny, nz = 3, 3, 3  # Number of unit cells in each direction
lattice, lattice_connectivity = duplicate_unit_cell(unit_cell, unit_cell_connectivity, nx, ny, nz)

#Quick check aka debugging
'''
print("lattice", lattice)
print(len(lattice))
print(type(lattice))
'''

# Visualize the lattice (optional)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot atoms
ax.scatter(lattice[:, 0], lattice[:, 1], lattice[:, 2], color='blue')

# Plot connections
for i, j in lattice_connectivity:
    ax.plot(*zip(lattice[i], lattice[j]), color='red')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()




#print("lattice", lattice)
latice_df = pd.DataFrame(lattice, columns=['X', 'Y', 'Z'])
#print(latice_df)
latice_df.to_csv('test_zebi_lattice.csv', index=False, header=False)
lattice_connectivity_df = pd.DataFrame(lattice_connectivity)
lattice_connectivity_df.to_csv('test_zebi_connectivity.csv', index=False, header=False)

list_length = len(lattice_connectivity)
print(list_length)
rayon = [0.5 for _ in range(2*list_length)]
rayon_df = pd.DataFrame(rayon)
print(rayon_df)
rayon_df.to_csv('test_zebi_rayon.csv', index=False, header=False)



