import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

'''Script error, the lattice is not generated correctly, the nodes are not connected correctly'''

#files
nodes_file = pd.read_csv('coordonnees_cube.csv')
nodes_file = np.concatenate((np.array([[0, 0, 0]]),nodes_file))
beams_file = pd.read_csv('connectivite_cube.csv') 
radius_file = pd.read_csv(r'CC\table_rayon_motif_CC.csv')


# Parameters
a = 1  # side length in mm

# Define the nodes and beams in the unit cell

test = np.array(nodes_file)

cell_nodes = np.array(nodes_file)
beam_cell = np.array(beams_file)

# Lattice dimensions
in_x, in_y, in_z = 2, 2, 2

# Generate the lattice
nodes = np.empty((0, 3))
beams = np.empty((0, 2), int)
for i in range(in_x):
    for j in range(in_y):
        for k in range(in_z):
            new_nodes = cell_nodes + np.array([i * a, j * a, k * a])
            new_beams = beam_cell + len(nodes)
            nodes = np.vstack((nodes, new_nodes))
            beams = np.vstack((beams, new_beams))

# Apply rotation (if necessary)
omega = 0  # rotation angle around the X-axis
M = np.array([
    [1, 0, 0],
    [0, np.cos(np.radians(omega)), -np.sin(np.radians(omega))],
    [0, np.sin(np.radians(omega)), np.cos(np.radians(omega))]
])
nodes = nodes @ M.T  # Apply rotation

# Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot beams
i = 0
for beam in beams:
    node1, node2 = nodes[beam[0]], nodes[beam[1]]
    ax.plot([node1[0], node2[0]], [node1[1], node2[1]], [node1[2], node2[2]], 'k')

# Plot nodes
ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], color='r')

# Set plot limits and labels
ax.set_xlim([0, a * in_x])
ax.set_ylim([0, a * in_y])
ax.set_zlim([0, a * in_z])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
