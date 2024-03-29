import sys
import numpy as np
import networkx as nx
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation

    
BOUNDING_BOX = [100] * 3


N_ASSOCIATIONS = 100
INLIERS = 1000
OUTLIERS = 10

MAX_TRANSLATION_AMOUNT = [100, 100, 100]


original = np.round(np.random.rand(N_ASSOCIATIONS, 3) * BOUNDING_BOX)

selected_indices = np.random.permutation(N_ASSOCIATIONS)[:INLIERS]
selected = np.vstack([
    original[selected_indices],
    np.random.rand(OUTLIERS, 3) * BOUNDING_BOX
])


t = np.random.rand(3) * MAX_TRANSLATION_AMOUNT
r = Rotation.from_euler('xyz', np.random.rand(3) * np.pi / 2)
R = r.as_matrix()


transformed = selected @ R.T + t

J = np.random.permutation(N_ASSOCIATIONS)
I = np.hstack([selected_indices, np.random.randint(0, N_ASSOCIATIONS, size=OUTLIERS)])[J]

original_distance_matrix = np.linalg.norm(original[I, None] - original[I], axis=2)
transformed_distance_matrix = np.linalg.norm(transformed[J, None] - transformed[J], axis=2)


m = 1
threshold = 0.1
consistency_matrix = np.exp(m * -np.abs(original_distance_matrix - transformed_distance_matrix))
consistency_matrix[consistency_matrix < threshold] = 0

X = np.arange(N_ASSOCIATIONS)
consistency_matrix[X, X] = 0

G = nx.from_numpy_array(consistency_matrix)
node_color = list(nx.core_number(G).values())

pose = nx.circular_layout(G)
nx.draw_networkx_nodes(G, pose, node_size=100, node_color=node_color)


edge_labels = nx.get_edge_attributes(G, "weight")
weight = list(edge_labels.values())

nx.draw_networkx_edges(G, pose, edgelist=edge_labels, width=1, alpha=weight)



plt.title(f'Density: {nx.density(G):0.2f}')
plt.show()


if len(sys.argv) > 1:
    with open(sys.argv[1], 'w') as f:
        f.write(f'{N_ASSOCIATIONS}\n')
        np.savetxt(f, consistency_matrix)
        J_inv = np.zeros_like(J)
        J_inv[J] = X
        true_association = J_inv[:INLIERS]
        f.write(f'{INLIERS}\n')
        np.savetxt(f, true_association)

