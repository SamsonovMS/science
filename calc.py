"""
Particle tracing.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode
from scipy.stats import maxwell


def init_xyz(n_part, box=[[15.0, 16.0], [-1.0, 1.0], [-1.0, 1.0]]):
    """
    Initial positions of particles.

    :param n_part: number of particles.
    :return: array of positions (n_part, 3).
    """
    ret = np.zeros((n_part, 3))
    thickness = [el[1] - el[0] for el in box]
    for i in range(3):
        for j, el in enumerate(np.random.rand(n_part)):
            ret[j][i] = el*thickness[i] + box[i][0]
    return ret


def init_vel(n_part, v_mean=[1, 0.062, 0.062]):
    """
    Maxwell for radial and perpendicular (thermal) initial velocities.
    Since OX points to the Sun, vx should be negative.
    Since mean value = 2*scale*sqrt(2/pi), I'll put v_mean/1.59577 as scale.

    :param n_part: number of particles.
    :param v_mean: list of 3 mean velocities x, y, z. 0.062 = 35/(400*sqrt(2))
    :return: array of velocities (n_part, 3).
    """
    ret = np.zeros((n_part, 3))
    for i in range(3):
        for j, el in enumerate(maxwell.rvs(scale=v_mean[i]/1.59577, size=n_part)):
            ret[j][i] = el
            if i == 0:
                ret[j][i] *= -1
            elif np.random.randint(2): # some particles should have negative vy or vz or both
                ret[j][i] *= -1
    return ret


def fun(t, _y, f_args):
    """
    Right hand side of the differential equations

    dx/dt = w_x
    dy/dt = w_y
    dz/dt = w_z
    dw_x/dt = (E_x + delta*(b_z*w_y - w_z*b_y))/theta
    dw_y/dt = (E_y + delta*(b_x*w_z - w_x*b_z))/theta
    dw_z/dt = (E_z + delta*(b_y*w_x - w_y*b_x))/theta
    
    Current sheet is rotated around the point (x_axis, z_axis) on the corner angle.
    """
    x, y, z, w_x, w_y, w_z = _y
    delta, theta, dh, dz, angle, x_axis, z_axis, case = f_args
    e_x = 0.0
    e_y = 0.0
    e_z = 0.0
    # take rotation into account
    r_point = np.sqrt(np.power(x - x_axis, 2) + np.power(z - z_axis, 2))
    a_point = np.arctan(float(z - z_axis) / (x - x_axis))
    z_cs = z_axis + r_point * np.sin(a_point - angle)
    b_x = 0.0
    b_y = 0.0
    b_z = np.tanh(x / dh)  # shock wave
    b_x += - case * np.tanh(z_cs / dz) # rotated current sheet
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


def traces(n_steps, n_part, t_koeff, *f_args):
    """
    Main function to calc trajectories.

    :param n_steps: number of time steps
    :param n_part: number of particles
    :param t_koeff: max time = n_steps*t_koeff
    :return: trajectories in 6D for all time steps
    """
    # Create an `ode` instance to solve the system of differential
    # equations defined by `fun`, and set the solver method.
    solver = ode(fun)
    solver.set_integrator('dopri5')
    # Give the value of delta and theta to the solver. This is passed to
    # `fun` when the solver calls it.
    solver.set_f_params(f_args)

    # Create the array `t` of time values at which to compute
    # the solution.
    max_time = n_steps * t_koeff
    t = np.linspace(0.0, max_time, n_steps)

    # Create an array to hold the solutions.
    sol = [np.empty((n_steps, 6)) for j in range(n_part)]

    # Initialize positions
    xyz = init_xyz(n_part)

    # Initialize maxwell velocities
    v = init_vel(n_part)

    for i in range(n_part):
        # Set the initial value _y(0) = _y0.
        _y0 = [[] for j in range(6)]
        _y0[:3] = xyz[i]
        _y0[3:] = v[i]
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

    plt.savefig('trajectories.png')
    plt.show()

    return sol


def show_hist(values, coord=0):
    """
    Show histogram of values.

    :param values: array(n_part, 3)
    """
    fig, ax = plt.subplots(1, 1)
    buf = []
    for el in values:
        buf.append(el[coord])
    ax.hist(buf, normed=False, histtype='stepfilled', alpha=0.2, bins=100)
    plt.show()


if __name__ == "__main__":
    show_hist(init_xyz(10000), 2)