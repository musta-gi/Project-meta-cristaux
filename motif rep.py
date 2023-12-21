import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Définition du motif initial CFC
points = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5), (1, 1, 1)]
connectivite = [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (3, 4)]

# Fonction pour répéter le motif
def repeter_motif(points, connectivite, repetitions, decalage):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for i in range(repetitions):
        for j in range(repetitions):
            for k in range(repetitions):
                # Décalage du motif
                decale_points = [(x + i * decalage[0], y + j * decalage[1], z + k * decalage[2]) for x, y, z in points]
                for p1, p2 in connectivite:
                    ax.plot([decale_points[p1][0], decale_points[p2][0]], 
                            [decale_points[p1][1], decale_points[p2][1]],
                            [decale_points[p1][2], decale_points[p2][2]], 'k-')

# Répéter le motif
repetitions = 2  # Nombre de fois à répéter le motif dans chaque direction
decalage = (1, 1, 1)  # Décalage à appliquer à chaque répétition
repeter_motif(points, connectivite, repetitions, decalage)

plt.show()
