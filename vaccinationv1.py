# Vaccination v1
# imports
from matplotlib import pyplot as plt
from vaccination import Vaccination
import numpy as np
from random import uniform
from math import ceil

# margin is a hyperparameter to be determined with the experiments
pi_margin = 0.05
delta_margin = 0.0125
low_center = 2*delta_margin
high_center = -low_center


############# Vaccination Rate Set Partitioning #######################
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

############# Vaccination Rate Change Set Partitioning ################


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

############################# HELPERS #################################


def maximum_membership(vaccination_rate):
    # Returns the class of the current vaccination rate status as an index
    memberships = [low_vac_rate(vaccination_rate), opt_vac_rate(
        vaccination_rate), high_vac_rate(vaccination_rate)]
    idx = memberships.index(max(memberships))
    membership = memberships[idx]
    return idx, membership


def MaxCriterion(index, membership):
    # Throws a random Vaccination rate change depending on the current
    #  vaccinating rate membership value to the corresponding set
    if index == 0:  # Low Vac Rate
        left = (membership-1)*delta_margin + low_center
        right = low_center - (membership-1) * delta_margin
        # left = delta_margin*(2*membership-1) + low_center
        # right = 0.25
    elif index == 1:  # Opt Vac Rate
        left = (membership-1)*delta_margin + 0
        right = 0 - (membership-1) * delta_margin
    else:  # High Vac Rate
        left = (membership-1)*delta_margin + high_center
        right = high_center - (membership-1) * delta_margin
        # left = -0.25
        # right = high_center + delta_margin*(1-2*membership)
    return round(uniform(left, right), 2)


def round_up_to_nearest_ten(x):
    # Just to be used in .viewVaccination()
    return int(ceil(x/10))*10


def find_stabilized_point(l):
    # Finds convergence point from vaccinated_percentage_curve_ attribute
    error = 0.01
    left = 0.6-error
    right = 0.6+error
    counter = 0
    found = False
    for i in range(len(l)):
        if left < l[i] < right:
            if not found:
                first_point = i
            found = True
            counter += 1
            if counter >= 30:
                return first_point
        else:
            found = False
            counter = 0
    return False


############################# VACCINATION ###############################

def Vaccinationv1():
    vaccination = Vaccination()
    days = 500
    convergence_day = 100

    for day in range(days):
        # First get pi
        pi, _ = vaccination.checkVaccinationStatus()
        # calculate membership and return set and the membership
        idx, membership = maximum_membership(pi)
        # calculate output control variable from membership using MaxCriterion
        delta = MaxCriterion(idx, membership)
        # Update pi_dot
        vaccination.vaccinatePeople(delta)

    # Once the loop finishes, get convergence day
    convergence_day = find_stabilized_point(
        vaccination.vaccinated_percentage_curve_)

    return vaccination, convergence_day


convergence_day = 0
# Visualize results
while True:
    vaccination, convergence_day = Vaccinationv1()
    if not convergence_day == 0 and convergence_day < 100:
        break

vaccination.viewVaccination(convergence_day, np.sum(
    vaccination.vaccination_rate_curve_[1:convergence_day]), save_dir="FIGSv1\\denemeler")
