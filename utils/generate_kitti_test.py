
import numpy as np
import open3d as o3d
import networkx as nx 
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

from tqdm import tqdm
import pyvista as pv

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--source', type=str, required=True)
parser.add_argument('--target', type=str, required=True)
parser.add_argument('--sample', type=int, default=1000)
parser.add_argument('--seed', type=int, required=False)

    

args = parser.parse_args()
if args.seed is not None:
    np.random.seed(args.seed)


print("PROCESSING SOURCE")
source = o3d.geometry.PointCloud()
source_points = np.load(args.source)
source.colors = o3d.utility.Vector3dVector([[1, 0, 0] for i in range(len(source_points))])
source.points = o3d.utility.Vector3dVector(source_points)

source.estimate_normals(o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=30))
source_fpfh = o3d.pipelines.registration.compute_fpfh_feature(source, o3d.geometry.KDTreeSearchParamHybrid(radius=20, max_nn=100))


print("PROCESSING TARGET")
target = o3d.geometry.PointCloud()
target_points = np.load(args.target)
target.points = o3d.utility.Vector3dVector(target_points)
target.colors = o3d.utility.Vector3dVector([[0, 0, 1] for i in range(len(target_points))])

target.estimate_normals(o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=30))
target_fpfh = o3d.pipelines.registration.compute_fpfh_feature(target, o3d.geometry.KDTreeSearchParamHybrid(radius=20, max_nn=100))


print("FIND MATCHES BETWEEN FPFH DESCRIPTORS")
tree = cKDTree(target_fpfh.data.T)
dists_st, nn_inds = tree.query(source_fpfh.data.T, k=1, workers=-1)

corres_st_ind_t = nn_inds
corres_st_ind_s = np.arange(len(nn_inds))


print("MUTUAL FILTER")
tree = cKDTree(source_fpfh.data.T)
dists_ts, nn_inds = tree.query(target_fpfh.data.T, k=1, workers=-1)

corres_ts_ind_s = nn_inds
corres_ts_ind_t = np.arange(len(nn_inds))


mutual_filter = (corres_ts_ind_s[corres_st_ind_t] == corres_st_ind_s)
mutual_filter[:] = False
mutual_filter[np.argsort(dists_st)[:1000]] = True

source_indices = corres_st_ind_s[mutual_filter]
target_indices = corres_st_ind_t[mutual_filter]
dists = dists_st[mutual_filter]

source_corr = source_points[source_indices]
target_corr = target_points[target_indices]


N = len(source_indices)
print(f'FPFH generates {N} putative correspondences.')


points = np.concatenate((source_corr, target_corr), axis=0)
lines = []
for i in range(N):
    lines.append([i, i + N])

colors = [[0, 1, 0] for i in range(len(lines))]

line_set = o3d.geometry.LineSet(
    points=o3d.utility.Vector3dVector(points),
    lines=o3d.utility.Vector2iVector(lines),
)
line_set.colors = o3d.utility.Vector3dVector(colors)

o3d.visualization.draw_geometries([source, target, line_set])


print("CREATING CONSISTENCY MATRIX")
source_distance_matrix = np.linalg.norm(source_points[source_indices, None] - source_points[source_indices], axis=2)
target_distance_matrix = np.linalg.norm(target_points[target_indices, None] - target_points[target_indices], axis=2)

consistency_matrix = np.exp(-
    np.square(source_distance_matrix - target_distance_matrix)
    #(np.linalg.norm(source_distance_matrix, axis=1) + np.linalg.norm(target_distance_matrix, axis=1))
)

consistency_matrix[consistency_matrix < 0.7] = 0

for x in np.unique(source_indices):
    i, = np.where(source_indices == x)
    consistency_matrix[np.ix_(i, i)] = 0


for x in np.unique(target_indices):
    j, = np.where(target_indices == x)
    consistency_matrix[np.ix_(j, j)] = 0


X = np.arange(len(target_indices))
consistency_matrix[X, X] = 0


print("FINDING MAXIMUM CLIQUE")
G = nx.from_numpy_array(consistency_matrix)

print(f"GRAPH DENSITY: {nx.density(G)}")
biggest = []
for clique in nx.clique.find_cliques(G):
    if len(clique) > len(biggest):
        print(f"FOUND A CLIQUE OF SIZE {len(clique)}: {clique}")
        biggest = clique
    

if args.output:
    with open(args.output, 'w') as f:
        f.write(f'{len(consistency_matrix)}\n')
        np.savetxt(f, consistency_matrix)
        true_association = biggest#max(cliques, key=len)
        f.write(f'{len(true_association)}\n')
        np.savetxt(f, true_association)
        print(f"WRITTEN TO {args.output}")


source_indices = source_indices[biggest]
target_indices = target_indices[biggest]

correspondance = np.vstack([source_indices, target_indices])
correspondance = o3d.utility.Vector2iVector(correspondance.T)


T = o3d.pipelines.registration.TransformationEstimationPointToPoint(False).compute_transformation(
    source,
    target,
    correspondance
)

o3d.visualization.draw_geometries([source.transform(T), target])


# thresholds = np.linspace(1, 0, 300)

# lens = []
# biggest = []

# pbar = tqdm(thresholds)
# for t in pbar:
#     consistency_matrix[consistency_matrix > t] = 0
#     G = nx.from_numpy_array(consistency_matrix)
#     cliques = list(nx.clique.find_cliques(G))
#     lens.append(len(cliques))
#     biggest.append(len(max(cliques, key=len)))
#     pbar.set_description(f'Threshold: {t:02f} - ({lens[-1], biggest[-1]})')


# fig, ax = plt.subplots()
# axes = [ax, ax.twinx()]
# plt1, = axes[0].plot(thresholds, lens)
# plt2, = axes[1].plot(thresholds, biggest, color='red')
# plt.legend([plt1, plt2], ['N. Maximal Cliques', 'Maximum Clique'], loc='lower right')
# plt.show()
