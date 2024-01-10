import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os

#---------------------------------------------
def extract_data(coordinates, connectivity):
    """Extract data from files."""
    cells_table = pd.read_csv(coordinates, header=None)
    beams_table = pd.read_csv(connectivity, header=None)
    cells = np.array(cells_table)
    beams = [tuple(x) for x in beams_table.values]
    return cells, beams


def generate_unit_cell(cells, beams, constant):
    """Generate a unit cell."""
    
    # Normalize the cells' values between 0 and 1
    normalized_cells = cells / np.max(cells)
    # Multiply by the constant
    atoms = normalized_cells * constant
    connectivity = beams
    return atoms, connectivity

def duplicate_unit_cell(unit_cell, connectivity, constant, nx, ny, nz):
    """
    Duplicate the unit cell to form a lattice with connectivity.

    Parameters:
    unit_cell (numpy.ndarray): Coordinates of atoms in the unit cell.
    connectivity (list): List of tuples representing connected atom indices.
    constant (float): Lattice constant, the length of the unit cell cube.
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
                translation = np.array([x, y, z]) * constant
                # Translate and add the unit cell to the lattice
                translated_cell = unit_cell + translation
                lattice.append(translated_cell)

                # Add connectivity for this unit cell
                offset = num_atoms_per_cell * (x * ny * nz + y * nz + z)
                for conn in connectivity:
                    lattice_connectivity.append((conn[0] + offset, conn[1] + offset))
    
    return np.vstack(lattice), lattice_connectivity

def save_to_DF(coordinates, connectivity, width):
    """Save the lattice to a DataFrame."""
    
    coordinates_DF = pd.DataFrame(coordinates, columns=['X', 'Y', 'Z'])
    connectivity_DF = pd.DataFrame(connectivity)
    radius_list_length = max(len(coordinates), len(connectivity))
    radius = [width for _ in range(radius_list_length)]
    radius_DF = pd.DataFrame(radius)
    
    return coordinates_DF, connectivity_DF, radius_DF
    

def plot_lattice(lattice, lattice_connectivity):
    """Plot the lattice."""
    
    # Visualize the lattice (optional)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot atoms
    #ax.scatter(lattice[:, 0], lattice[:, 1], lattice[:, 2], color='blue')
    translated_lattice, translated_connectivity = translate_lattice(lattice, lattice_connectivity)
    # Plot connections
    for i, j in lattice_connectivity:
        ax.plot(*zip(lattice[i], lattice[j]), color='red')
    

    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

def export_lattice(name, destination, coordinates_DF, connectivity_DF, radius_DF):
    """Export the lattice to a csv file."""
    
    folder_path = os.path.join(destination, name)
    os.makedirs(folder_path, exist_ok=True)
    
    coordinates_DF.to_csv(os.path.join(folder_path, 'coordinates.csv'), index=False, header=False)
    connectivity_DF.to_csv(os.path.join(folder_path, 'connectivity.csv'), index=False, header=False)
    radius_DF.to_csv(os.path.join(folder_path, 'radius.csv'), index=False, header=False)
#---------------------------------------------
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def rotate_unit_cell(lattice, connectivity, axis, angle):
    """
    Rotate the lattice about the given axis by the given angle.

    Parameters:
    lattice (numpy.ndarray): Coordinates of atoms in the lattice.
    connectivity (list): List of tuples representing connected atom indices.
    axis (numpy.ndarray): Axis of rotation.
    angle (float): Angle of rotation in radians.

    Returns:
    numpy.ndarray: Coordinates of atoms in the rotated lattice.
    list: Connectivity of the rotated lattice.
    """
    
    # Rotation matrix
    rotation = rotation_matrix(axis, angle)

    
    rotated_cell = []
    for point in lattice:
        rotated_cell.append(np.dot(point, rotation))
    rotated_cell = np.array(rotated_cell)
    
    '''# translate cell to remove negative values
    translated_cell = []
    for point in rotated_cell:
        translated_cell.append(point + np.abs(np.min(rotated_cell)))
    # normalize cell
    rotated_cell = rotated_cell / np.max(rotated_cell)
    #print("rotated_cell", rotated_cell)'''
    
    return rotated_cell, connectivity

def translate_lattice(lattice, lattice_connectivity):
    coord_min = np.min(lattice)
    coord_max = np.max(lattice)
    translation_constant = np.sqrt(2)*12 #np.sqrt(2)*(coord_max - coord_min)
    
    print("translation_constant", translation_constant)
    
    translated_cell = []
    for point in lattice:
        translated_cell.append(point + translation_constant)
    
    translated_cell = np.array(translated_cell)
    return translated_cell, lattice_connectivity

def plan_equation(point1, point2, point3):
    """Return the equation of the plane defined by 3 points."""
    
    #print("point1", point1)
    #print("point2", point2)
    #print("point3", point3)
    
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    x3, y3, z3 = point3
    
    a = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1)
    b = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1)
    c = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
    d = -(a*x1 + b*y1 + c*z1)
    
    return a, b, c, d

def which_side(point, plane):
    """Return the side of the plane on which the point is."""
    
    a, b, c, d = plane
    x, y, z = point
    
    result = a*x + b*y + c*z + d
    
    if result > 0:
        return 1
    elif result < 0:
        return -1
    else:
        return 0
    
def cut_lattice(lattice, lattice_connectivity, plane, side):
    """Cut the lattice with the given plane.""" # adjust this with a hash table
    cut_lattice = []
    cut_connectivity = []
    
    for i, point in enumerate(lattice):
        if which_side(point, plane) == side:
            cut_lattice.append(point)
    
    for i, connection in enumerate(lattice_connectivity):
        if which_side(lattice[connection[0]], plane) == side and which_side(lattice[connection[1]], plane) == side:
            cut_connectivity.append(connection)
    
    return np.array(cut_lattice), cut_connectivity
    

def cut_lattice():
    pass
    
#---------------------------------------------
if __name__ =='__main__':
    
    coordinates = 'CFC\coordonnees_motif_CFC.csv'
    connectivity = 'CFC\connectivite_motif_CFC.csv'
    constant = 6
    width = 0.5
    n = 2
    nx, ny, nz = n, n, n  # Number of unit cells in each direction
    name = 'test_G7R'
    destination = ''
    
    axis = [1, 1, 1]
    angle = np.pi/3
    
    cells, beams = extract_data(coordinates, connectivity)
    atoms, connectivity = generate_unit_cell(cells, beams, constant)
    
    
    lattice, lattice_connectivity = duplicate_unit_cell(atoms, connectivity, constant, nx, ny, nz)
    
    atoms, connectivity = rotate_unit_cell(lattice, lattice_connectivity, axis, angle) #rotate the generated lattice
    coordinates_DF, connectivity_DF, radius_DF = save_to_DF(atoms, connectivity, width)
    plot_lattice(atoms, connectivity)
    #plot_lattice(lattice, lattice_connectivity)
    export_lattice(name, destination, coordinates_DF, connectivity_DF, radius_DF)