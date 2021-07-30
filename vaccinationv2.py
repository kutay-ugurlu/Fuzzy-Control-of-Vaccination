# Vaccination v2
from matplotlib import pyplot as plt
from vaccination import Vaccination
from random import uniform
import numpy as np
from math import ceil

# margin is a hyperparameter to be determined with the experiments
pi_margin = 0.05
pi_dot_margin = 0.1
delta_margin = 0.025
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

####### Input Vaccination Rate Change Set Partitioning #######################


def low_vac_rate_input(x):
    if x < 0 - pi_dot_margin:
        return 1
    elif x > 0:
        return 0
    else:
        return -1/pi_dot_margin * x + 0/pi_dot_margin


def high_vac_rate_input(x):
    if x < 0:
        return 0
    elif x > 0 + pi_dot_margin:
        return 1
    else:
        return 1/pi_dot_margin * x - 0/pi_dot_margin


def opt_vac_rate_input(x):
    if x < 0-pi_dot_margin or x > 0+pi_dot_margin:
        return 0
    else:
        return np.sign(x-0) * ((0-x)/pi_dot_margin) + 1

############# Output Vaccination Rate Change Set Partitioning #######################


def increase_low(x, low_center):
    if x < low_center-delta_margin or x > low_center+delta_margin:
        return 0
    else:
        return np.sign(x-low_center) * ((low_center-x)/delta_margin) + 1


def decrease_low(x, high_center):
    if x < high_center-delta_margin or x > high_center+delta_margin:
        return 0
    else:
        return np.sign(x-high_center) * ((high_center-x)/delta_margin) + 1


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


def opt_vac_rate_change(x):
    if x < 0-delta_margin or x > 0+delta_margin:
        return 0
    else:
        return np.sign(x-0) * ((0-x)/delta_margin) + 1

############################# HELPERS #################################


def maximum_membership_pi(vaccination_rate):
    # Returns the class of the current vaccination rate status as an index
    memberships = [low_vac_rate(vaccination_rate), opt_vac_rate(
        vaccination_rate), high_vac_rate(vaccination_rate)]
    idx = memberships.index(max(memberships))
    membership = memberships[idx]
    return idx, membership


def maximum_membership_pi_dot(vaccination_rate):
    # Returns the class of the current vaccination rate change status as an index
    memberships = [low_vac_rate_input(vaccination_rate), opt_vac_rate_input(
        vaccination_rate), high_vac_rate_input(vaccination_rate)]
    idx = memberships.index(max(memberships))
    membership = memberships[idx]
    return idx, membership


def MaxCriterion(pi_idx, pi_dot_idx, membership):
    # Throws a random Vaccination rate per day change depending on the current
    #  vaccinating rate membership value to the corresponding set
    if pi_idx == 0:  # Low Vac Rate
        if pi_dot_idx == 0 or pi_dot_idx == 1:  # Low OR Opt
            # INCREASE HIGH
            left = delta_margin*(2*membership-1) + low_center
            right = 0.25
        else:
            # INCREASE LOW
            left = (membership-1)*delta_margin + low_center
            right = low_center - (membership-1) * delta_margin

    elif pi_idx == 2:  # High Vac Rate
        if pi_dot_idx == 1 or pi_dot_idx == 2:  # High OR Opt
            # DECREASE HIGH
            left = -0.25
            right = high_center + delta_margin*(1-2*membership)
        else:
            # DECREASE LOW
            left = (membership-1)*delta_margin + low_center
            right = low_center - (membership-1) * delta_margin
    else:  # Opt Vac Rate
        if pi_dot_idx == 0 or pi_dot_idx == 1:  # Low OR Opt
            # OPTIMAL CHANGE
            left = (membership-1)*delta_margin
            right = -(membership-1) * delta_margin
        else:
            # DECREASE HIGH
            left = -0.25
            right = high_center + delta_margin*(1-2*membership)

    # Return random number from alpha cut
    return round(uniform(left, right), 3)


def round_up_to_nearest_ten(x):
    return int(ceil(x/10))*10


def find_stabilized_point(l):
    # Finds convergence point from vaccinated_percentage_curve_ attribute
    error = 0.005
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
            if counter >= 60:
                return first_point
        else:
            found = False
            counter = 0
    return False


############################# VACCINATION ###############################

def vaccinationv2():
    vaccination = Vaccination()

    days = 200
    convergence_day = 200

    for day in range(days):
        # First get pi and pi_dot
        current, effective = vaccination.checkVaccinationStatus()
        # calculate memberships and return sets and the memberships
        pi_dot_idx, _ = maximum_membership_pi_dot(effective)
        pi_idx, membership = maximum_membership_pi(current)
        # calculate output control variable from pi and pi_dot indices and pi membership using MaxCriterion
        delta = MaxCriterion(pi_idx, pi_dot_idx, membership)
        # Update pi_dot
        vaccination.vaccinatePeople(delta)

    # Once the loop finishes, get convergence day
    convergence_day = find_stabilized_point(
        vaccination.vaccinated_percentage_curve_)
    return convergence_day, vaccination


convergence_day = 0
# Visualize results
while True:
    convergence_day, vaccination = vaccinationv2()
    if not convergence_day == 0 and convergence_day < 20:
        break

vaccination.viewVaccination(convergence_day, np.sum(
    vaccination.vaccination_rate_curve_[1:convergence_day]), save_dir="FIGSv2\\denemeler\\")
