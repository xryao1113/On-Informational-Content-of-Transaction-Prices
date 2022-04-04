import numpy as np
import matplotlib.pyplot as plt


def Exogenous_PPC(S, n, c):
    S1_payoffs = []
    S2_payoffs = []
    for i in range(S):
        buyer_types = list(np.random.uniform(0, 1, 1)) * n
        theta_M = max(buyer_types)
        if c <= 0.25:
            S1_payoff = theta_M - c
            S2_payoff = theta_M
        else:
            if theta_M >= 0.5:
                S1_payoff = 0.5
                S2_payoff = 0.5
            else:
                S1_payoff = 0
                if theta_M >= 0.25:
                    S2_payoff = 0.25
                else:
                    S2_payoff = 0
        S1_payoffs.append(S1_payoff)
        S2_payoffs.append(S2_payoff)
    return S1_payoffs, S2_payoffs


def Exogenous_ID(S, n, c):
    c1 = n / (n + 1) * (1 - (1 / (n + 1)) ** (1 / n))

    S1_payoffs = []
    S2_payoffs = []
    for i in range(S):
        buyer_types = list(np.random.uniform(0, 1, n))
        buyer_types.sort()
        theta_M = buyer_types[-1]
        theta_m = buyer_types[-2]
        if c <= c1:
            S1_payoff = theta_M - c
            c2 = (n - 1) * theta_M / n * (1 - n ** (1 / (1 - n)))
            p2 = theta_M * (n ** (1 / (1 - n)))
            if c <= c2:
                S2_payoff = theta_m - c
            elif theta_m >= p2:
                S2_payoff = p2
            else:
                S2_payoff = 0
        else:
            p1 = (1 / (n + 1)) ** (1 / n)
            if theta_M >= p1:
                S1_payoff = p1
                p2 = ((n + 1) ** (1 / n) - 1) * (n * (n + 1)) ** (1 /
                                                                  (1 - n)) + (2 - (n + 1) ** (1 / n)) * n ** (1 / (1 - n))
                if theta_m >= p2:
                    S2_payoff = p2
                else:
                    S2_payoff = 0
            else:
                S1_payoff = 0
                p2 = (1 / (n + 1)) ** (2 / n)
                if theta_M >= p2:
                    S2_payoff = p2
                else:
                    S2_payoff = 0
        S1_payoffs.append(S1_payoff)
        S2_payoffs.append(S2_payoff)
    return S1_payoffs, S2_payoffs


def Inifinite_Horizon(S, n, c, delta):
    S1_payoffs = []
    S2_payoffs = []

    for i in range(S):
        buyer_types = list(np.random.uniform(0, 1, 1)) * n
        theta_M = max(buyer_types)

        set_probability = (1 - delta) * (1 - 2 * c) / (2 * c * delta)
        set_probability = min(set_probability, 1)
        set_probability = max(set_probability, 0)

        stage_counter = 0
        while True:
            if np.random.uniform(0, 1, 1) < set_probability:
                S1_payoffs.append((theta_M - c) * delta ** stage_counter)
                S2_payoffs.append(theta_M * delta ** stage_counter)
                break
            elif np.random.uniform(0, 1, 1) < set_probability:
                S1_payoffs.append(theta_M * delta ** stage_counter)
                S2_payoffs.append((theta_M - c) * delta ** stage_counter)
                break
            stage_counter += 1
    return S1_payoffs, S2_payoffs


np.random.seed(0)
S = 5000
n = 100
c_range = np.arange(0.01, 1, 0.01)

plt.figure(figsize=(15, 10))
S1_payoff_means = []
S2_payoff_means = []
for c in c_range:
    S1_payoffs, S2_payoffs = Exogenous_PPC(S, n, c)
    S1_payoff_means.append(np.mean(S1_payoffs))
    S2_payoff_means.append(np.mean(S2_payoffs))
plt.plot(c_range, S1_payoff_means, c='crimson', label='S1')
plt.plot(c_range, S2_payoff_means, c='royalblue', label='S2')
plt.xlabel('Cost of Market Research')
plt.ylabel('Average Utility')
plt.xlim([0.01, 1])
plt.ylim([0, 1])
plt.legend()
plt.savefig('Exogenous_PPC.png',
            facecolor="white", bbox_inches='tight')

plt.figure(figsize=(15, 10))
S1_payoff_means = []
S2_payoff_means = []
for c in c_range:
    S1_payoffs, S2_payoffs = Exogenous_ID(S, n, c)
    S1_payoff_means.append(np.mean(S1_payoffs))
    S2_payoff_means.append(np.mean(S2_payoffs))
plt.plot(c_range, S1_payoff_means, c='crimson', label='S1')
plt.plot(c_range, S2_payoff_means, c='royalblue', label='S2')
plt.xlabel('Cost of Market Research')
plt.ylabel('Average Utility')
plt.xlim([0.01, 1])
plt.ylim([0, 1])
plt.legend()
plt.savefig('Exogenous_ID.png',
            facecolor="white", bbox_inches='tight')

c_range = np.arange(0.01, 0.25, 0.01)
delta_range = [0.2, 0.5, 0.8]
plt.figure(figsize=(15, 10))
for delta in delta_range:
    S1_payoff_means = []
    S2_payoff_means = []
    for c in c_range:
        S1_payoffs, S2_payoffs = Inifinite_Horizon(S, n, c, delta)
        S1_payoff_means.append(np.mean(S1_payoffs))
        S2_payoff_means.append(np.mean(S2_payoffs))
    plt.plot(c_range, S1_payoff_means, marker='s',
             label='S1, delta = ' + str(delta))
    plt.plot(c_range, S2_payoff_means, marker='o',
             label='S2, delta = ' + str(delta))
plt.xlabel('Cost of Market Research')
plt.ylabel('Average Utility')
plt.xlim([0.01, 0.25])
plt.ylim([0, 1])
plt.legend()
plt.savefig('Infinite_Horizon.png',
            facecolor="white", bbox_inches='tight')
