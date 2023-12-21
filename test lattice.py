import numpy as np
import pandas as pd

# Parameters
a = 10  # side length in mm

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

# Writing to a file
with open('lattice_structure.txt', 'w') as file:
    file.write('*node\n')
    for i, node in enumerate(nodes, start=1):
        file.write(f'{node[0]}, {node[1]}, {node[2]}\n')

    file.write('*element\n')
    for i, beam in enumerate(beams, start=1):
        file.write(f'{beam[0] + 1}, {beam[1] + 1}\n')

nodes_df = pd.DataFrame(nodes)
nodes_df.to_csv('nodes.csv', index=False)
beams_df = pd.DataFrame(beams)
beams_df.to_csv('beams.csv', index=False)
rayon = [0.5 for i in range(len(beams))]
rayon_df = pd.DataFrame(rayon)
rayon_df.to_csv('rayon.csv', index=False)
print("File 'lattice_structure.txt' has been created.")
