"""
Particle tracing.
See full explanation https://docs.google.com/document/d/
1Uebu5muzgRV0MCNWihRWR2x5Dqxgk6k5Sj-HnUkoc_M/edit#heading=h.51eqdmgucqiq
"""

import numpy as np
import matplotlib.pyplot as plt
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
    """
    x, y, z, w_x, w_y, w_z = _y
    delta, theta, dh, d, dx1, dx2, bz0, XO, XM, XNL, bz2 = f_args
    e_x = 0.0
    e_y = 0.0
    e_z = 0.0
    b_x = np.tanh(z / d)
    b_y = 0.0
    b_z = np.tanh(x / dh)
    if x < XM:
        b_z += -bz0 * np.tanh((x - XO) / dx1)
    else:
        b_z += bz2 * np.tanh((x - XNL) / dx2)
    return [w_x, w_y, w_z,
            (e_x + delta * (b_z * w_y - w_z * b_y)) / theta,
            (e_y + delta * (b_x * w_z - w_x * b_z)) / theta,
            (e_z + delta * (b_y * w_x - w_y * b_x)) / theta]


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
bz2 = bz0 * np.tanh((XM - XO) / dx1) / np.tanh((XNL - XM) / dx2)
f_args = (delta, theta, dh, d, dx1, dx2, bz0, XO, XM, XNL, bz2)
solver.set_f_params(f_args)

# Number of time steps.
n_steps = 300

# Create the array `t` of time values at which to compute
# the solution.
t0 = 0.0
t1 = 500
t = np.linspace(t0, t1, n_steps)

# Number of particles.
n_part = 10

# Create an array to hold the solutions.
sol = [np.empty((n_steps, 6)) for j in range(n_part)]

for i in range(n_part):
    # Set the initial value _y(0) = _y0.
    _y0 = [[] for j in range(6)]
    _y0[:3] = 10 * (np.random.random_sample(3) - np.random.random_sample(3))
    _y0[3:] = (np.random.random_sample(3) - np.random.random_sample(3)) / 10
    solver.set_initial_value(_y0, t0)

    # Put the initial value in the solutions array.
    sol[i][0] = _y0

    # Repeatedly call the `integrate` method to advance the
    # solution to time t[k], and save the solution in sol[k].
    k = 1
    while solver.successful() and solver.t < t1:
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
