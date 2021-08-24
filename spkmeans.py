import sys
import numpy as np
import pandas as pd
import myspkmeans


def pad(x):
    if x >= 0:
        return str(x) + (6 - len(str(x))) * "0"
    else:
        return str(x) + (7 - len(str(x))) * "0"


try:
    k = int(sys.argv[1])
    if k < 0:
        print("Invalid Input!")
        exit()
except ValueError:
    print("Invalid Input!")
    exit()
except IndexError:
    print("Invalid Input!")
    exit()

allowed_goals = ["spk", "wam", "ddg", "lnorm", "jacobi"]

try:
    goal = str(sys.argv[2])
    if goal not in allowed_goals:
        print("Invalid Input!")
        exit()
except IndexError:
    print("Invalid Input!")
    exit()

try:
    filename = str(sys.argv[3])
except IndexError:
    print("Invalid Input!")
    exit()

data = myspkmeans.fit(filename, goal, k)
if goal != "spk":
    exit()
# print(data)
# print("----")
df = pd.DataFrame(data)
# print(df)
# print("---")
dim = len(df.columns)
# print("dim = ",dim)
# print("----")
k = dim
df['cluster'] = -1.0
points_to_cluster_init = df.to_numpy().flatten().tolist()
del df['cluster']

np_df = df.to_numpy()
num_of_points = len(np_df)

df['d'] = np.nan
df['p'] = np.nan
np_df_ext = df.to_numpy()

np.random.seed(0)
rand_start = np.random.randint(0, (len(df) - 1))
centroids_df_indices = [rand_start]
z = 1

while z < k:
    for i in range(num_of_points):
        ind = centroids_df_indices[z - 1]
        d_cand = (np.linalg.norm(np_df[i] - np_df[ind])) ** 2
        if z == 1:
            np_df_ext[i][dim] = d_cand
        else:
            if d_cand < np_df_ext[i][dim]:
                np_df_ext[i][dim] = d_cand
    d_sum = np_df_ext.sum(axis=0)[dim]
    for i in range(num_of_points):
        np_df_ext[i][dim + 1] = np_df_ext[i][dim] / d_sum

    n_array = np.arange(num_of_points)
    centroids_df_indices.append(np.random.choice(n_array, p=np_df_ext[:, dim + 1]))
    z += 1

centroids_lists = [np_df[ind].tolist() for ind in centroids_df_indices]  # need to output this to file
centroids_init = [num for sublist in centroids_lists for num in sublist]

to_print = ','.join([str(num) for num in centroids_df_indices])
print(to_print)

centroids = myspkmeans.fit2(k, num_of_points, dim, centroids_init, points_to_cluster_init,
                            len(centroids_init), len(points_to_cluster_init))

centroids = np.asarray(centroids).round(4)
centroids = np.array_split(centroids, k)
centroids = [centroids[i].tolist() for i in range(len(centroids))]

for cent in centroids:
    print(','.join(pad(val) for val in cent))
