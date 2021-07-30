from matplotlib import pyplot as plt
from vaccination import Vaccination
import numpy as np

# margin is a hyperparameter to be determined with the experiments
pi_margin = 0.05
delta_margin = 0.0125
low_center = 2*delta_margin
high_center = -low_center


####### Vaccination Rate Set Partitioning #######################
def low_vac_rate(x):
    if x < 0.6 - pi_margin:
        return 1
    elif x > 0.6:
        return 0
    else:
        return -1/pi_margin * x + 0.6/pi_margin


def high_vac_rate(x):
    if x < 0.6:
        return 0
    elif x > 0.6 + pi_margin:
        return 1
    else:
        return 1/pi_margin * x - 0.6/pi_margin


def opt_vac_rate(x):
    if x < 0.6-pi_margin or x > 0.6+pi_margin:
        return 0
    else:
        return np.sign(x-0.6) * ((0.6-x)/pi_margin) + 1

############# Vaccination Rate Change Set Partitioning #######################


def low_vac_rate_change(x, low_center):
    if x < low_center-delta_margin or x > low_center+delta_margin:
        return 0
    else:
        return np.sign(x-low_center) * ((low_center-x)/delta_margin) + 1


def high_vac_rate_change(x, high_center):
    if x < high_center-delta_margin or x > high_center+delta_margin:
        return 0
    else:
        return np.sign(x-high_center) * ((high_center-x)/delta_margin) + 1


def opt_vac_rate_change(x):
    if x < 0-delta_margin or x > 0+delta_margin:
        return 0
    else:
        return np.sign(x-0) * ((0-x)/delta_margin) + 1


def increase_high(x, high_center):
    if x > high_center+delta_margin:
        return 1
    elif x < high_center - delta_margin:
        return 0
    else:
        return 1/(2*delta_margin) * x - (high_center-delta_margin)/(2*delta_margin)


def decrease_high(x, low_center):
    if x < low_center - delta_margin:
        return 1
    elif x > low_center + delta_margin:
        return 0
    else:
        return -1/(2*delta_margin) * x + (low_center+delta_margin)/(2*delta_margin)

##########################################################################


plt.subplot(1, 2, 1)
pi = np.arange(0, 1, 0.05)
low_rate = list(map(low_vac_rate, pi))
high_rate = list(map(high_vac_rate, pi))
opt_rate = list(map(opt_vac_rate, pi))
plt.plot(pi, low_rate)
plt.plot(pi, opt_rate)
plt.plot(pi, high_rate)
plt.title("Fuzzy Partition of the Sets for margin="+str(pi_margin))
plt.xlabel("$\pi$")
plt.ylabel("Membership Values")
plt.legend(["Low Vaccination", "Optimal Vaccination", "High Vaccination"],)


plt.subplot(1, 2, 2)
delta = np.arange(-0.25, 0.25, 0.01)
low_rate = list(map(lambda item: increase_high(item, low_center), delta))
high_rate = list(
    map(lambda item: decrease_high(item, high_center), delta))
opt_rate = list(map(opt_vac_rate_change, delta))
plt.plot(delta, low_rate)
plt.plot(delta, opt_rate)
plt.plot(delta, high_rate)
plt.title("Fuzzy Partition of the Sets for margin="+str(delta_margin))
plt.xlabel("$\delta$")
plt.ylabel("Membership Values")
plt.legend(["Low Vaccination \n Rate Change",
            "Optimal Vaccination \n Rate Change", "High Vaccination \n Rate Change"],)
plt.show()
