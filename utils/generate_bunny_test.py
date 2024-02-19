import os
import numpy as np
import open3d as o3d
import networkx as nx
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation

import urllib.request
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--rho', type=float, required=True)
parser.add_argument('--assoc', type=int, required=True)
parser.add_argument('--output', action='store_true', required=False)
parser.add_argument('--outliers', type=int, default=0, required=False)
parser.add_argument('--noise', type=float, default=0.5, required=False)

args = parser.parse_args()
N_ASSOCIATIONS = args.assoc
OUTLIER_RATIO = float(args.rho)
NOISE = args.noise

N_OUTLIERS = round(OUTLIER_RATIO * N_ASSOCIATIONS)
N_INLIERS = N_ASSOCIATIONS - N_OUTLIERS


FILE = "bunny.pcd"
if not os.path.exists(FILE):
    URL = "https://raw.githubusercontent.com/PointCloudLibrary/pcl/master/test/bunny.pcd"
    contents = urllib.request.urlretrieve(URL, FILE)

original = np.asarray(o3d.io.read_point_cloud(FILE).points) * 1000
original = original[np.random.permutation(len(original))]

R = Rotation.random().as_matrix()
transformed = original @ R.T + np.max(original) + np.random.rand(*original.shape) * NOISE

if args.outliers is not None:
    original = np.vstack([
        original,
        np.random.rand(args.outliers, 3) * 1000
    ])


I = np.hstack([np.arange(N_INLIERS), np.random.randint(N_INLIERS, len(original), size=N_OUTLIERS)])
J = np.hstack([np.arange(N_INLIERS), np.random.randint(N_INLIERS, len(transformed), size=N_OUTLIERS)])

P = np.random.permutation(N_ASSOCIATIONS)
I = I[P]
J = J[P]


original_distance_matrix = np.linalg.norm(original[I, None] - original[I], axis=2)
transformed_distance_matrix = np.linalg.norm(transformed[J, None] - transformed[J], axis=2)

threshold = 0.1
consistency_matrix = np.exp(-np.abs(original_distance_matrix - transformed_distance_matrix))
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

if args.output:
    with open(f'bunny-{args.rho * 100:0.0f}-{args.assoc}.txt', 'w') as f:
        f.write(f'{N_ASSOCIATIONS}\n')
        np.savetxt(f, consistency_matrix)
        true_association = np.where(P < N_INLIERS)[0]
        f.write(f'{N_INLIERS}\n')
        np.savetxt(f, true_association)
        print(f"FILE WIRTTEN TO bunny-{args.rho * 100:0.0f}-{args.assoc}.txt")