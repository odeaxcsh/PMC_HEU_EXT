import os
import numpy as np
import open3d as o3d
import networkx as nx
import matplotlib.pylab as plt
from scipy.spatial.transform import Rotation

from tqdm import tqdm
import urllib.request
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--rho', type=float, required=True)
parser.add_argument('--assoc', type=int, required=True)
parser.add_argument('--output', action='store_true', required=False)
parser.add_argument('--outliers', type=int, default=0, required=False)
parser.add_argument('--noise', type=float, default=0, required=False)
parser.add_argument('--threshold', type=float, default=0.5, required=False)
parser.add_argument('-q', '--quite', action='store_true', required=False)
parser.add_argument('--seed', type=int, required=False)


def expected_distance_to_noise(sigma, bins=100, n=1000000):
    sample = np.random.normal(scale=sigma * np.sqrt(1/2), size=n)
    sample = np.abs(sample)
    values, base = np.histogram(sample, bins=bins)
    base = np.convolve(base, np.ones(2) / 2, 'valid')
    values = values / n
    cumulative = np.cumsum(values)
    return base, cumulative
    

args = parser.parse_args()
if args.seed is not None:
    np.random.seed(args.seed)


N_ASSOCIATIONS = args.assoc
OUTLIER_RATIO = args.rho
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


I = np.hstack([np.arange(N_INLIERS), np.random.randint(len(original), size=N_OUTLIERS)])
J = np.hstack([np.arange(N_INLIERS), np.random.randint(len(transformed), size=N_OUTLIERS)])

Z, = np.where((I == J))
Z = Z[Z >= N_INLIERS]
I[Z] += 1


P = np.random.permutation(N_ASSOCIATIONS)
I = I[P]
J = J[P]


original_distance_matrix = np.linalg.norm(original[I, None] - original[I], axis=2)
transformed_distance_matrix = np.linalg.norm(transformed[J, None] - transformed[J], axis=2)

consistency_matrix = np.exp(-np.square(original_distance_matrix - transformed_distance_matrix))


# clean consistency matrix:
for x in np.unique(I):
    i, = np.where(I == x)
    consistency_matrix[np.ix_(i, i)] = 0


for x in np.unique(J):
    j, = np.where(J == x)
    consistency_matrix[np.ix_(j, j)] = 0


X = np.arange(N_ASSOCIATIONS)
consistency_matrix[X, X] = 0
print("CONSISTENCY MATRIX CREATED")

consistency_matrix[consistency_matrix < args.threshold] = 0

G = nx.from_numpy_array(consistency_matrix)

if not args.quite:
    node_color = list(nx.core_number(G).values())

    pose = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pose, node_size=100, node_color=node_color, cmap='Reds')


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


cliques = list(nx.clique.find_cliques(G))
print(f"N. Maximal Cliques: {len(cliques)}")

    # limit = 1
    # while len(cliques) > 1000:
    #     limit += 1
    #     cliques = list(filter(lambda a: len(a) > limit, cliques))
    # cliques = list(sorted(cliques, key=len, reverse=True))[:1000]

    # print(f"N. Maximal Cliques (|C| > {limit}): {len(cliques)}")

    # shared = np.zeros((len(cliques), len(cliques)))

    # for i, c in tqdm(enumerate(cliques)):
    #     for j, cp in enumerate(cliques):
    #         if i == j:
    #             continue
    #         shared[i, j] = len(set(c) & set(cp))


    # plt.imshow(shared)
    # plt.colorbar()
    # plt.show()
#     maximums.append(len(max(cliques, key=len)))
#     shareds.append(shared / (maximums[-1] - 1))


# fig, axes = plt.subplots(4, 5)
# axes = axes.ravel()

# for i in range(20):
#     img = axes[i].imshow(shareds[i], vmin=0, vmax=max(s.max() for s in shareds))
#     axes[i].set_title(f"rho: {rhos[i]:0.02f}")
#     axes[i].set_xticks([], [])
#     axes[i].set_yticks([], [])
    
# fig.colorbar(img, ax=axes.tolist())
# plt.show()
# consistency_matrix = consistency_matrix_org

# noise_distro = expected_distance_to_noise(NOISE)
# thresholds = np.linspace(noise_distro[0].max(), 0, 100)

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
# axes = [ax, ax.twinx(), ax.twinx()]
# plt1, = axes[0].plot(thresholds, lens)
# plt2, = axes[1].plot(thresholds, biggest, color='red')
# axes[2].set_yticks([], [])
# axes[2].plot(*noise_distro, alpha=0.1)
# plt.legend([plt1, plt2], ['N. Maximal Cliques', 'Maximum Clique'], loc='lower right')
# plt.show()
