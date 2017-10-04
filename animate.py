import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import pandas as pd
from os import chdir, makedirs
import datetime
# my modules
import calc
import convert
import fields
import c
import info


# Make new folder for each run
name = str(datetime.datetime.now()).split('.')[0]
makedirs(c.path + name)
chdir(c.path + name)

# Constants
case = 1
num = [400, 0, 200]
min_x = -30
max_x = 60
min_z = -25
max_z = 40
x = np.linspace(min_x, max_x, num=num[0])
z = np.linspace(min_z, max_z, num=num[2])

n_steps = 10
n_part = 2 # > 1
t_koeff = 0.5 # bigger => more calculations
delta = 1.0
theta = 1.0

dh = 0.75
d = 0.75
dx1 = 5.0
dx2 = 1.0
bz0 = 1.0
xo = 7.0
xm = 12.0
xnl = 15.0
bz2 = bz0 * np.tanh((xm - xo) / dx1) / np.tanh((xnl - xm) / dx2)

angle = np.pi/4
x_axis = -1.0
z_axis = 0.0

# Write constants to constants.txt
info.write_constants(case=case, num=num, min_x=min_x, max_x=max_x, min_z=min_z, max_z=max_z,
                     n_steps=n_steps, n_part=n_part, t_koeff=t_koeff, delta=delta, theta=theta,
                     dh=dh, d=d, dx1=dx1, dx2=dx2, bz0=bz0, xo=xo, xm=xm, xnl=xnl, bz2=bz2,
                     angle=angle, x_axis=x_axis, z_axis=z_axis)

# Calc, show and save png of total field of shock wave and rotated plasmoid
u, v, m = fields.main(num, x, z, case, dh, d, dx1, dx2, bz0, xo, xm, xnl, bz2, angle, x_axis, z_axis)
fields.quiver_(x, z, u, v, m)

# Calc, show and save png of trajectories
sol = calc.traces(case, n_steps, n_part, t_koeff, delta, theta, dh, d, dx1, dx2, bz0, xo, xm, xnl, bz2,
                  angle, x_axis, z_axis, case)
xyz = convert.sol_to_xyz(sol, n_steps, n_part)

# Animation
t = np.array([np.ones(n_part)*i for i in range(n_steps)]).flatten()
df = pd.DataFrame({'time': t , 'x': xyz[:, 0], 'y': xyz[:, 1], 'z': xyz[:, 2]})


def update_graph(num):
    data = df[df['time'] == num]
    graph.set_data(data.x, data.y)
    graph.set_3d_properties(data.z)
    title.set_text('step={}'.format(num))
    return title, graph


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('')

data = df[df['time'] == 0]
graph, = ax.plot(data.x, data.y, data.z, linestyle="", marker="o")

anim = matplotlib.animation.FuncAnimation(fig, update_graph, n_steps,
                               interval=1, blit=False)

xx, yy, zz = convert.field_to_surface(num, x, z, m)
ax.contourf(xx, yy, zz, zdir='y', offset=0, alpha=0.25, cmap=cm.coolwarm)

ax.view_init(elev=0., azim=-90)
ax.set_xlabel('X')
ax.set_xlim(-20, 20)
ax.set_ylabel('Y')
ax.set_ylim(-10, 10)
ax.set_zlabel('Z')
ax.set_zlim(-10, 10)

plt.show()