# Import Python modules
import numpy as np  # Version 1.20.3
import matplotlib.pyplot as plt  # Version 3.4.3


def Threshold_Difference(n):
    """
    Compute difference between two sides of the inequality (refer to Theorem 5.4)

    - n: the number of buyers
    """
    # Threshold for S2 (left-hand side of the inequality)
    S_2_threshold = (n - 1) / n * (1 - n ** (1 / (1 - n)) * (n + 2
                                                             - (1 / (n + 1)) ** (1 / n) - (n ** 2 - 1) / (n + 1 - (n + 1) ** (1 / n))))
    # Threshold for S1 (right-hand side of the inequality)
    S_1_threshold = n / (n + 1) * (1 - (1 / (n + 1)) ** (1 / n))

    return S_1_threshold - S_2_threshold  # Return the difference


# Specify range of n to compute difference on
n_range = list(np.arange(2, 1 * 10 ** 4, 1))
# Iterate over n_range to compute and store differences
difference_list = [Threshold_Difference(i) for i in n_range]

print('All difference between thresholds are greater than zero:',
      all(difference > 0 for difference in difference_list))

# Plot figure and save as png file in current directory
plt.figure(figsize=(30, 10), dpi=100)
plt.plot(n_range, difference_list)
plt.savefig('difference.png',
            facecolor="white", bbox_inches='tight')
