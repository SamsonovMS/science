"""
Particle tracing.
See
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode


def fun(t, _y, f_args):
    """
    Right hand side of the differential equations

    dx/dt = w_x
    dy/dt = w_y
    dz/dt = w_z
    dw_x/dt = (E_x + delta*(b_z*w_y - w_z*b_y))/theta
    dw_y/dt = (E_y + delta*(b_x*w_z - w_x*b_z))/theta
    dw_z/dt = (E_z + delta*(b_y*w_x - w_y*b_x))/theta
    
    Plasmoid is rotated around a point (x_axis, z_axis) on the corner angle.
    """
    x, y, z, w_x, w_y, w_z = _y
    delta, theta, dh, d, dx1, dx2, bz0, XO, XM, XNL, bz2, angle, x_axis, z_axis, case = f_args
    e_x = 0.0
    e_y = 0.0
    e_z = 0.0
    # take rotation into account
    r_point = np.sqrt(np.power(x - x_axis, 2) + np.power(z - z_axis, 2))
    a_point = np.arctan(float(z - z_axis) / (x - x_axis))
    x_plasm = x_axis + r_point * np.cos(a_point - angle)
    z_plasm = z_axis + r_point * np.sin(a_point - angle)
    b_x = case * np.tanh(float(z_plasm) / d)
    b_y = 0.0
    b_z = np.tanh(x / dh)  # shock wave
    if x_plasm < XM: # change +/- here
        b_z += - case * bz0 * np.tanh((x_plasm - XO) / dx1)
    else:
        b_z += case * bz2 * np.tanh((x_plasm - XNL) / dx2)
    r_value = np.sqrt(np.power(b_x, 2) + np.power(b_z, 2))
    a_value = np.arctan(b_z / b_x)
    if b_x < 0:
        a_value = np.pi + a_value
    b_x = r_value*np.cos(a_value + angle)
    b_z = r_value*np.sin(a_value + angle)
    return [w_x, w_y, w_z,
            (e_x + delta * (b_z * w_y - w_z * b_y)) / theta,
            (e_y + delta * (b_x * w_z - w_x * b_z)) / theta,
            (e_z + delta * (b_y * w_x - w_y * b_x)) / theta]


def traces(case, n_steps=300, n_part=1, t_coeff=0.5):
    """
    Main function to calc trajectories.

    :param n_steps: number of time steps
    :param n_part: number of particles
    :param t_coeff: max time = n_steps*t_coeff
    :return: trajectories in 6D for all time steps
    """
    # Create an `ode` instance to solve the system of differential
    # equations defined by `fun`, and set the solver method.
    solver = ode(fun)
    solver.set_integrator('dopri5')
    # Give the value of delta and theta to the solver. This is passed to
    # `fun` when the solver calls it.
    delta = 1.0
    theta = 1.0
    dh = 0.75
    d = 0.75
    dx1 = 5.0
    dx2 = 1.0
    bz0 = 1.0  # ?
    XO = 7.0  # ?
    XM = 12.0
    XNL = 15.0
    angle = 0.24
    x_axis = -6.0
    z_axis = 0.0
    bz2 = bz0 * np.tanh((XM - XO) / dx1) / np.tanh((XNL - XM) / dx2)
    f_args = (delta, theta, dh, d, dx1, dx2, bz0, XO, XM, XNL, bz2, angle, x_axis, z_axis, case)
    solver.set_f_params(f_args)

    # Create the array `t` of time values at which to compute
    # the solution.
    max_time = n_steps*t_coeff
    t = np.linspace(0.0, max_time, n_steps)

    # Create an array to hold the solutions.
    sol = [np.empty((n_steps, 6)) for j in range(n_part)]

    for i in range(n_part):
        # Set the initial value _y(0) = _y0.
        _y0 = [[] for j in range(6)]
        _y0[:3] = 10 * (np.random.random_sample(3) - np.random.random_sample(3))
        _y0[3] = (np.random.random_sample(1)*4 - np.random.random_sample(1)) / 10
        _y0[4:] = (np.random.random_sample(2) - np.random.random_sample(2)) / 10
        solver.set_initial_value(_y0, 0.0)

        # Put the initial value in the solutions array.
        sol[i][0] = _y0

        # Repeatedly call the `integrate` method to advance the
        # solution to time t[k], and save the solution in sol[k].
        k = 1
        while solver.successful() and solver.t < max_time:
            solver.integrate(t[k])
            sol[i][k] = solver.y
            k += 1


    # Plot the solutions.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(n_part):
        ax.scatter(sol[i][:, 0][::2], sol[i][:, 1][::2], sol[i][:, 2][::2], s=0.05)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

    return sol