import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# Example data (replace with your 5x5 results)
# rows = initial conditions
# cols = repeated runs
# -------------------------------------------------

# data = np.array([
#     [902621, 910083, 922876, 903261, 924899],
#     [912003, 892173, 909313, 896631, 922907],
#     [923071, 910079, 923793, 914806, 912328],
#     [923988, 899351, 902004, 917330, 921761],
#     [924196, 911640, 918768, 893970, 906958]
# ])

data = np.array([
    [3138431, 3120265, 3128086, 3150801, 3087835],
    [3148623, 3182599, 3221874, 3200711, 3094341],
    [3159512, 3193476, 3127823, 3115213, 3023976],
    [3070809, 3076668, 3158798, 3137750, 3015883],
    [3104320, 3183080, 3173591, 3166955, 3109303]
])

n_ic, n_runs = data.shape

# -------------------------------------------------
# Mean calculations
# -------------------------------------------------

grand_mean = np.mean(data)
condition_means = np.mean(data, axis=1)

# -------------------------------------------------
# Within-condition variance (GPU stochasticity)
# -------------------------------------------------

within_var = np.mean(np.var(data, axis=1, ddof=1))

# -------------------------------------------------
# Between-condition variance
# -------------------------------------------------

between_var = np.var(condition_means, ddof=1)

# -------------------------------------------------
# Derived metrics
# -------------------------------------------------

within_std = np.sqrt(within_var)
between_std = np.sqrt(between_var)

std_ratio = within_std / between_std

icc = between_var / (between_var + within_var)

# -------------------------------------------------
# Print results
# -------------------------------------------------

print("Grand mean:", grand_mean)

print("\nBetween-condition variance (initial condition effect):")
print("Variance:", between_var)
print("Std:", between_std)

print("\nWithin-condition variance (GPU stochasticity):")
print("Variance:", within_var)
print("Std:", within_std)

print("\nStd ratio (run / initial condition):", std_ratio)

cv_run = within_std / grand_mean
cv_ic  = between_std / grand_mean

print("CV run:", cv_run)
print("CV initial condition:", cv_ic)

print("\nIntraclass correlation (ICC):", icc)

# -------------------------------------------------
# Visualization
# -------------------------------------------------

plt.figure()
plt.boxplot(data.T)
plt.xlabel("Initial Condition")
plt.ylabel("Result Value")
plt.title("Run Variability per Initial Condition")
plt.show()

plt.figure()
plt.bar(["Initial Condition", "Run Stochasticity"],
        [between_std, within_std])
plt.ylabel("Standard Deviation")
plt.title("Standard Deviation Decomposition")
plt.show()