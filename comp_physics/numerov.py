import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import brentq

def f(current, E):
    return -2*(1/current+E)

def w_aux(current, h, E, f):
    return 1 - h**2*f(current+h, E)/12

def dx_aux(current, h, E, f):
    return 1 - h**2*f(current+h, E)/6

def xh_aux(current, h, E, f):
    return 2+5*h**2*f(current, E)/6

def numerator(current, h, x_current, dx_current, E, f):
    return xh_aux(current, h, E, f) * w_aux(current, -h, E, f) * x_current + 2*h*dx_current*dx_aux(current, -h, E, f)

def denominator(current, h, x_current, dx_current, E, f):
    return w_aux(current, h, E, f) * dx_aux(current, -h, E, f) + w_aux(current, -h, E, f) * dx_aux(current, h, E, f)

def xh(current, h, x_current, dx_current, E, f):
    return numerator(current, h, x_current, dx_current, E, f) / denominator(current, h, x_current, dx_current, E, f)


def numerov_log(E, f, current, start=1e-20, end=200, x_current=0, dx_current=0.0007, direction=-1, n_steps=5000):
    # E = -0.1
    # current = 10
    # h = 0.00001
    # x_current = 0
    # dx_current = 0.0007

    # x_prev = xh(current, -h, x_current, dx_current)


    if direction not in [-1, 1]:
        raise Exception("Invalid direction")

    # n_steps = int(current / h)

    grids = np.logspace(np.log(start), np.log(current), n_steps+1, base=np.e)[:]
    # grids = np.linspace(1e-5, current, n_steps+2)[1:]

    r_values = np.zeros((n_steps, ), dtype=float)
    xh_values = np.zeros((n_steps, ), dtype=float)


    # for i in tqdm(range(n_steps), position=0, leave=True):
    for i in range(n_steps):

        h = grids[::direction][i+1] -  grids[::direction][i]
        # print("h: ", h)

        current = grids[::direction][i]

        x_next = xh(current, h, x_current, dx_current, E, f)

        x_next_leap_frog = xh(current, 0.75*h, x_current, dx_current, E, f)

        dx_current = (x_next - x_next_leap_frog)/(0.25*h)

        x_current = x_next
        # current = current + h * direction

        r_values[::direction][i] = current
        xh_values[::direction][i] = x_current

    normalization_factor = np.sqrt((xh_values * xh_values).sum())

    xh_values /= normalization_factor

    return r_values, xh_values*direction, xh_values[0], xh_values[-1]


def find_eigen(f, E_lower=-0.7, E_upper=-1e-5, n_levels=100, current=200, start=1e-20, end=200, x_current=0, dx_current=0.0007, direction=-1, n_steps=5000):

    if E_lower * E_upper > 0:
        energy_grids = np.sign(E_lower) * np.logspace(np.log(abs(E_lower)), np.log(abs(E_upper)), n_levels, base=np.e)
    else:
        energy_grids = np.linspace(E_lower, E_upper, n_levels)

    E_guessed = []
    xh_guessed = []
    E_last = E_lower

    find_root = lambda e: numerov_log(e, f, current, start=start, x_current=x_current, dx_current=dx_current, direction=-1, n_steps=n_steps)[2]
    _, _, u0_last, _ = numerov_log(E_last, f,  current, start=start, x_current=x_current, dx_current=dx_current, direction=-1, n_steps=n_steps)
    # for egs in tqdm(E_guess, position=0, leave=True):
    for egs in energy_grids:
        rv, xhv, u0_current, _ = numerov_log(egs, f, current, start=start, x_current=x_current, dx_current=dx_current, direction=-1, n_steps=n_steps)
        # print("egs: ", egs, "u_0: ", u0_current)

        if u0_current * u0_last < 0:
            print("guess range: ", E_last, egs)

            egs_ = brentq(find_root, E_last, egs)


            print("energy level: ", egs_)

            E_guessed.append(egs_)
            xh_guessed.append((rv, xhv))

        E_last = egs

        u0_last = u0_current

    return E_guessed, xh_guessed


if __name__ == "__main__":
    current = 200
    # h = 0.00001
    x_current = 0
    dx_current = 0.0007
    start = 1e-20
    n_steps = 5000

    find_eigen(f, start=start, n_steps=n_steps)



