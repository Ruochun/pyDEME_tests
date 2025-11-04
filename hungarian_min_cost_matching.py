import numpy as np
from scipy.optimize import linear_sum_assignment

import matplotlib.pyplot as plt

# number of repeated slumping runs
num_runs = 10
# number of permutated runs
num_permutations = 10

diff_runs = []
total_dist = []
# compared all runs wrt to No.1
A = np.loadtxt(f"slump_va_final_1.csv", delimiter=',', comments='#')
total_dist.append(np.sum(np.linalg.norm(A, axis=1)))
for i in range(2, num_runs+1):
    B = np.loadtxt(f"slump_va_final_{i}.csv", delimiter=',', comments='#')

    cost = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    row_ind, col_ind = linear_sum_assignment(cost)
    total_distance = cost[row_ind, col_ind].sum()
    diff_runs.append(total_distance)
    dist = 0.
    for p in B:
        dist += np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
    total_dist.append(dist)

    # print(f"The mapping is {list(zip(row_ind, col_ind))}")

print(f"The total distances among runs are {diff_runs}")
print(f"The total distances of all particles is {total_dist}")

diff_permu = []
# compared all permutated runs wrt to No.1
A = np.loadtxt(f"slump_va_final_1.csv", delimiter=',', comments='#')
for i in range(2, num_permutations+1):
    B = np.loadtxt(f"slump_va_final_permu_{i}.csv", delimiter=',', comments='#')

    cost = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    row_ind, col_ind = linear_sum_assignment(cost)
    total_distance = cost[row_ind, col_ind].sum()
    diff_permu.append(total_distance)

print(f"The total distances among permutated runs are {diff_permu}")

# plot both in the same figure as segmented lines to show comparison
plt.figure(figsize=(8, 6))
plt.plot(range(1, num_runs), diff_runs, marker='o', label='Repeated Runs')
plt.plot(range(1, num_permutations), diff_permu, marker='s', label='Permutated Runs')
plt.xlabel('Run Index (compared to Run 1)')
plt.ylabel('Total Distance')
plt.title('Comparison of Total Distances Among Runs')
plt.legend()
# y limit is 0 to max + 10%
max_val = max(max(diff_runs), max(diff_permu))
plt.ylim(0, max_val * 1.1)
plt.grid(True)
plt.savefig('comparison_of_runs.png', dpi=300)
plt.show()
