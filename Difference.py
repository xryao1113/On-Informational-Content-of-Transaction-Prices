import numpy as np
import matplotlib.pyplot as plt


def threshold_difference(n):
    S_2_threshold = (n - 1) / n * (1 - n ** (1 / (1 - n)) * (n + 2
                                                             - (1 / (n + 1)) ** (1 / n) - (n ** 2 - 1) / (n + 1 - (n + 1) ** (1 / n))))
    S_1_threshold = n / (n + 1) * (1 - (1 / (n + 1)) ** (1 / n))

    return S_1_threshold - S_2_threshold


n_range = list(np.arange(2, 1 * 10 ** 4, 1))
difference_list = [threshold_difference(i) for i in n_range]

print('All difference between thresholds are greater than zero:',
      all(difference > 0 for difference in difference_list))

plt.figure(figsize=(30, 10), dpi=100)
plt.plot(n_range, difference_list)
