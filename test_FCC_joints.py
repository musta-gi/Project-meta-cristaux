import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

def generate_fcc_unit_cell(a):
    """
    Generate a face-centered cubic (FCC) unit cell.

    Parameters:
    a (float): Lattice constant, the length of the unit cell cube.

    Returns:
    numpy.ndarray: Coordinates of atoms in the unit cell.
    """
    # Coordinates for the corners and face centers of the cube
    corners = a * np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                            [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])
    face_centers = a * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                                 [0.5, 0.5, 1], [0.5, 1, 0.5], [1, 0.5, 0.5]])
    
    atoms = np.vstack((corners, face_centers))

    # Define connectivity within the unit cell
    connectivity = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if np.linalg.norm(atoms[i] - atoms[j]) < 1.1 * a / np.sqrt(2):  # threshold for nearest neighbors
                connectivity.append((i, j))
    
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
    
def rotate_unit_cell(unit_cell, connectivity, axis, angle):
    """
    Rotate the unit cell about the given axis by the given angle.

    Parameters:
    unit_cell (numpy.ndarray): Coordinates of atoms in the unit cell.
    connectivity (list): List of tuples representing connected atom indices.
    axis (numpy.ndarray): Axis of rotation.
    angle (float): Angle of rotation in radians.

    Returns:
    numpy.ndarray: Coordinates of atoms in the rotated unit cell.
    list: Connectivity of the rotated unit cell.
    """
    # Rotation matrix
    rotation = np.array([[np.cos(angle) + axis[0] ** 2 * (1 - np.cos(angle)),
                          axis[0] * axis[1] * (1 - np.cos(angle)) - axis[2] * np.sin(angle),
                          axis[0] * axis[2] * (1 - np.cos(angle)) + axis[1] * np.sin(angle)],
                         [axis[1] * axis[0] * (1 - np.cos(angle)) + axis[2] * np.sin(angle),
                          np.cos(angle) + axis[1] ** 2 * (1 - np.cos(angle)),
                          axis[1] * axis[2] * (1 - np.cos(angle)) - axis[0] * np.sin(angle)],
                         [axis[2] * axis[0] * (1 - np.cos(angle)) - axis[1] * np.sin(angle),
                          axis[2] * axis[1] * (1 - np.cos(angle)) + axis[0] * np.sin(angle),
                          np.cos(angle) + axis[2] ** 2 * (1 - np.cos(angle))]])

    # Rotate the unit cell
    rotated_cell = np.dot(unit_cell, rotation)
    #print("rotated_cell", rotated_cell)

    # Rotate the connectivity
    rotated_connectivity = []
    for conn in connectivity:
        rotated_connectivity.append((np.dot(unit_cell[conn[0]], rotation), np.dot(unit_cell[conn[1]], rotation)))
    #print("rotated_connectivity", rotated_connectivity)
    
    return rotated_cell, connectivity

# Lattice constant (you can adjust this)
a = 20

# Generate the unit cell and its connectivity
unit_cell, unit_cell_connectivity = generate_fcc_unit_cell(a)

# Rotate the unit cell (optional)
unit_cell, unit_cell_connectivity = rotate_unit_cell(unit_cell, unit_cell_connectivity, axis=[1, 1, 1], angle=np.pi / 4)

# Duplicate the unit cell to create a lattice (change nx, ny, nz as needed)
nx, ny, nz = 2,2,2  # Number of unit cells in each direction
lattice, lattice_connectivity = duplicate_unit_cell(unit_cell, unit_cell_connectivity, nx, ny, nz)

print("lattice", lattice)
print("lattice_connectivity", lattice_connectivity)

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
latice_df.to_csv('joints_test_zebi_lattice.csv', index=False, header=False)
lattice_connectivity_df = pd.DataFrame(lattice_connectivity)
lattice_connectivity_df.to_csv('joints_test_zebi_connectivity.csv', index=False, header=False)

list_length = len(lattice_connectivity)
print(list_length)
rayon = [0.5 for _ in range(list_length)]
rayon_df = pd.DataFrame(rayon)
print(rayon_df)
rayon_df.to_csv('joints_test_zebi_rayon.csv', index=False, header=False)
