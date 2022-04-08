# Import Python modules
import numpy as np  # Version 1.20.3
import matplotlib.pyplot as plt  # Version 3.4.3


def Exogenous_PPC(S, n, c):
    """
    Numerical simulation of the informational model with exogenous seller orderings
    and buyer types of perfectly positive correlation.

    - S: the number of simulations
    - n: the number of buyers
    - c: the cost of market research
    """
    S1_payoffs = []
    S2_payoffs = []
    # Loop over simulations
    for i in range(S):
        # Random generation of buyer types
        buyer_types = list(np.random.uniform(0, 1, 1)) * n
        theta_M = max(buyer_types)

        # Case when S1 conducts research
        if c <= 0.25:
            S1_payoff = theta_M - c
            S2_payoff = theta_M
        # Case when S1 does not conduct research
        else:
            # Case when first transaction is successful
            if theta_M >= 0.5:
                S1_payoff = 0.5
                S2_payoff = 0.5
            # Case when first transaction is unsuccessful
            else:
                S1_payoff = 0
                # Case when second transaction is successful
                if theta_M >= 0.25:
                    S2_payoff = 0.25
                # Case when second transaction is unsuccessful
                else:
                    S2_payoff = 0

        # Record seller utilities
        S1_payoffs.append(S1_payoff)
        S2_payoffs.append(S2_payoff)

    return S1_payoffs, S2_payoffs


def Exogenous_ID(S, n, c):
    """
    Numerical simulation of the informational model with exogenous seller orderings
    and independently distributed buyer types.

    - S: the number of simulations
    - n: the number of buyers
    - c: the cost of market research
    """
    # Threshold of S1 for cost (refer to Lemma 5.3 for derivation)
    c1 = n / (n + 1) * (1 - (1 / (n + 1)) ** (1 / n))

    S1_payoffs = []
    S2_payoffs = []
    # Loop over simulations
    for i in range(S):
        # Random generation of buyer types
        buyer_types = list(np.random.uniform(0, 1, n))
        buyer_types.sort()
        theta_M = buyer_types[-1]
        theta_m = buyer_types[-2]

        # Case when S1 conducts research
        if c <= c1:
            S1_payoff = theta_M - c
            # Threshold of S2 for cost (refer to Theorem 5.4 for derivation)
            c2 = (n - 1) * theta_M / n * (1 - n ** (1 / (1 - n)))
            p2 = theta_M * (n ** (1 / (1 - n)))

            # Case when S2 conducts research
            if c <= c2:
                S2_payoff = theta_m - c
            # Case when second transaction is successful
            elif theta_m >= p2:
                S2_payoff = p2
            # Case when second transaction is unsuccessful
            else:
                S2_payoff = 0
        # Case when S1 does not conduct research
        else:
            p1 = (1 / (n + 1)) ** (1 / n)
            # Case when first transaction is successful
            if theta_M >= p1:
                S1_payoff = p1
                p2 = ((n + 1) ** (1 / n) - 1) * (n * (n + 1)) ** (1 /
                                                                  (1 - n)) + (2 - (n + 1) ** (1 / n)) * n ** (1 / (1 - n))
                # Case when second transaction is successful
                if theta_m >= p2:
                    S2_payoff = p2
                # Case when second transaction is unsuccessful
                else:
                    S2_payoff = 0
            # Case when first transaction is unsuccessful
            else:
                S1_payoff = 0
                p2 = (1 / (n + 1)) ** (2 / n)
                # Case when second transaction is successful
                if theta_M >= p2:
                    S2_payoff = p2
                # Case when second transaction is unsuccessful
                else:
                    S2_payoff = 0

        # Record seller utilities
        S1_payoffs.append(S1_payoff)
        S2_payoffs.append(S2_payoff)

    return S1_payoffs, S2_payoffs


def Inifinite_Horizon(S, n, c, delta):
    """
    Numerical simulation of the informational model with endogenous seller orderings
    and buyer types of perfectly positive correlation, along the infinite horizon.

    - S: the number of simulations
    - n: the number of buyers
    - c: the cost of market research
    - delta: the discount factor
    """
    S1_payoffs = []
    S2_payoffs = []
    # Loop over simulations
    for i in range(S):
        # Random generation of buyer types
        buyer_types = list(np.random.uniform(0, 1, 1)) * n
        theta_M = max(buyer_types)

        # first case of mixed strategies (refer to Theorem 8.1)
        if c <= 0.25:
            # Compute probability parameter
            set_probability = (1 - delta) * (1 - 2 * c) / (2 * c * delta)
            set_probability = min(set_probability, 1)
            set_probability = max(set_probability, 0)

            stage_counter = 0
            while True:  # Game continues as long as sellers are stalling
                # Case when S1 sets his price in the current stage
                if np.random.uniform(0, 1, 1) < set_probability:
                    S1_payoffs.append((theta_M - c) * delta ** stage_counter)
                    S2_payoffs.append(theta_M * delta ** stage_counter)
                    break
                # Case when S2 sets his price in the current stage
                elif np.random.uniform(0, 1, 1) < set_probability:
                    S1_payoffs.append(theta_M * delta ** (stage_counter + 1))
                    S2_payoffs.append((theta_M - c) * delta ** stage_counter)
                    break
                stage_counter += 1
        # second case of mixed strategies (refer to Theorem 8.1)
        else:
            # Compute probability parameter
            set_probability = 4 * (1 - delta) / delta
            set_probability = min(set_probability, 1)
            set_probability = max(set_probability, 0)

            stage_counter = 0
            while True:  # Game continues as long as sellers are stalling
                # Case when S1 sets his price in the current stage
                if np.random.uniform(0, 1, 1) < set_probability:
                    # Case when first transaction is successful
                    if theta_M >= 0.5:
                        S1_payoffs.append(0.5 * delta ** stage_counter)
                        S2_payoffs.append(0.5 * delta ** stage_counter)
                    # Case when first transaction is unsuccessful
                    else:
                        S1_payoffs.append(0)
                        # Case when second transaction is successful
                        if theta_M >= 0.25:
                            S2_payoffs.append(0.25 * delta ** stage_counter)
                        # Case when second transaction is unsuccessful
                        else:
                            S2_payoffs.append(0)
                    break
                # Case when S2 sets his price in the current stage
                elif np.random.uniform(0, 1, 1) < set_probability:
                    # Case when first transaction is successful
                    if theta_M >= 0.5:
                        S1_payoffs.append(0.5 * delta ** (stage_counter + 1))
                        S2_payoffs.append(0.5 * delta ** stage_counter)
                    # Case when first transaction is unsuccessful
                    else:
                        S2_payoffs.append(0)
                        # Case when second transaction is successful
                        if theta_M >= 0.25:
                            S1_payoffs.append(
                                0.25 * delta ** (stage_counter + 1))
                        # Case when second transaction is unsuccessful
                        else:
                            S1_payoffs.append(0)
                    break
                stage_counter += 1

    return S1_payoffs, S2_payoffs


# Set up of parametrs
np.random.seed(0)
S = 5000
n = 100
c_range = np.arange(0.01, 1, 0.01)
delta_range = [0.2, 0.5, 0.8]

# Plot of average utility from Exogenous_PPC over c_range
plt.figure(figsize=(15, 10))
S1_payoff_means = []
S2_payoff_means = []
# Loop over c values
for c in c_range:
    S1_payoffs, S2_payoffs = Exogenous_PPC(S, n, c)
    S1_payoff_means.append(np.mean(S1_payoffs))
    S2_payoff_means.append(np.mean(S2_payoffs))

# Plot configurations
plt.plot(c_range, S1_payoff_means, c='crimson', label='S1')
plt.plot(c_range, S2_payoff_means, c='royalblue', label='S2')
plt.xlabel('Cost of Market Research')
plt.ylabel('Average Utility')
plt.xlim([0.01, 1])
plt.ylim([0, 1])
plt.legend()
# Save plot as png file in current directory
plt.savefig('Exogenous_PPC.png',
            facecolor="white", bbox_inches='tight')

# Plot of average utility from Exogenous_ID over c_range
plt.figure(figsize=(15, 10))
S1_payoff_means = []
S2_payoff_means = []
# Loop over c values
for c in c_range:
    S1_payoffs, S2_payoffs = Exogenous_ID(S, n, c)
    S1_payoff_means.append(np.mean(S1_payoffs))
    S2_payoff_means.append(np.mean(S2_payoffs))

# Plot configurations
plt.plot(c_range, S1_payoff_means, c='crimson', label='S1')
plt.plot(c_range, S2_payoff_means, c='royalblue', label='S2')
plt.xlabel('Cost of Market Research')
plt.ylabel('Average Utility')
plt.xlim([0.01, 1])
plt.ylim([0, 1])
plt.legend()
# Save plot as png file in current directory
plt.savefig('Exogenous_ID.png',
            facecolor="white", bbox_inches='tight')

# Plot of average utilities from Infinite_Horizon over c_range and delta_range
plt.figure(figsize=(15, 10))
# Loop over delta values
for delta in delta_range:
    S1_payoff_means = []
    S2_payoff_means = []
    # Loop over c values
    for c in c_range:
        S1_payoffs, S2_payoffs = Inifinite_Horizon(S, n, c, delta)
        S1_payoff_means.append(np.mean(S1_payoffs))
        S2_payoff_means.append(np.mean(S2_payoffs))
    plt.plot(c_range, S1_payoff_means, marker='s',
             label='S1, delta = ' + str(delta))
    plt.plot(c_range, S2_payoff_means, marker='o',
             label='S2, delta = ' + str(delta))

# Plot configurations
plt.xlabel('Cost of Market Research')
plt.ylabel('Average Utility')
plt.xlim([0.01, 1])
plt.ylim([0, 1])
plt.legend()
# Save plot as png file in current directory
plt.savefig('Infinite_Horizon.png',
            facecolor="white", bbox_inches='tight')
