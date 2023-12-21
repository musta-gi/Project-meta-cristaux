import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

'''Translation of the Matlab script to Python, the lattice is not generated correctly, the nodes are not connected correctly'''

# Parameters
a = 1  # side length in mm

# Define the nodes and beams in the unit cell
cell_nodes = np.array([
    [0, 0, 0], [a, 0, 0], [a, a, 0], [0, a, 0],
    [a/2, a/2, 0], [0, a/2, a/2], [a/2, 0, a/2],
    [a/2, a, a/2], [a, a/2, a/2], [0, 0, a], [a, 0, a],
    [a, a, a], [0, a, a], [a/2, a/2, a]
])

beam_cell = np.array([
    [1, 5], [1, 6], [1, 7], [2, 5], [2, 7], [2, 9], [3, 5], [3, 8], [3, 9],
    [4, 5], [4, 6], [4, 8], [5, 6], [5, 7], [5, 8], [5, 9], [6, 7], [6, 8],
    [7, 9], [8, 9], [6, 14], [7, 14], [8, 14], [9, 14], [10, 14], [10, 6],
    [10, 7], [11, 14], [11, 7], [11, 9], [12, 14], [12, 8], [12, 9], [13, 14],
    [13, 6], [13, 8]
]) - 1  # Adjust for zero-indexing in Python

test_nodes_df = pd.DataFrame(cell_nodes)
test_nodes_df.to_csv('test_nodes.csv', index=False)
test_beams_df = pd.DataFrame(beam_cell)
test_beams_df.to_csv('test_beams.csv', index=False)

# Lattice dimensions
in_x, in_y, in_z = 1, 1, 1

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
