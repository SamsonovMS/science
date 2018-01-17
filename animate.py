# coding=utf-8
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
min_x = -10
max_x = 6
min_z = -4
max_z = 4
x = np.linspace(min_x, max_x, num=num[0])
z = np.linspace(min_z, max_z, num=num[2])

n_steps = 300
n_part = 50  # > 1
t_koeff = 0.5 # time of particle movement for one step. Bigger — larger trajectories and more rude results
delta = 0.245
theta = 1.0 # ion mass / proton mass
dh = 0.085
dz = 0.25
angle = -np.pi/6
x_axis = -1.0
z_axis = 0.0

# Write constants to constants.txt
info.write_constants(case=case, num=num, min_x=min_x, max_x=max_x, min_z=min_z, max_z=max_z,
                     n_steps=n_steps, n_part=n_part, t_koeff=t_koeff, delta=delta, theta=theta,
                     dh=dh, dz=dz, angle=angle, x_axis=x_axis, z_axis=z_axis)

# Calc, show and save png of total field of shock wave and rotated plasmoid
u, v, m = fields.main(x, z, case, dh, dz, angle, x_axis, z_axis)
fields.quiver_(x, z, u, v, m)

# Calc, show and save png of trajectories
sol = calc.traces(case, n_steps, n_part, t_koeff, delta, theta, dh, dz,
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