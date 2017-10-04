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

case = 1
num = [400, 0, 200]
x = np.linspace(-30, 60, num=num[0])
z = np.linspace(-25, 40, num=num[2])

u, v, m = fields.main(num, x, z, case)
fields.quiver_(x, z, u, v, m)

xx, yy, zz = convert.field_to_surface(num, x, z, m)

n_steps = 30
n_part = 5 # > 1

a = convert.sol_to_xyz(calc.traces(case, n_steps=n_steps, n_part=n_part), n_steps, n_part)

t = np.array([np.ones(n_part)*i for i in range(n_steps)]).flatten()
df = pd.DataFrame({'time': t , 'x': a[:, 0], 'y': a[:, 1], 'z': a[:, 2]})


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

ani = matplotlib.animation.FuncAnimation(fig, update_graph, n_steps - 1,
                               interval=1, blit=False)
ax.contourf(xx, yy, zz, zdir='y', offset=0, alpha=0.25, cmap=cm.coolwarm)

ax.view_init(elev=0., azim=-90)
ax.set_xlabel('X')
ax.set_xlim(-20, 20)
ax.set_ylabel('Y')
ax.set_ylim(-10, 10)
ax.set_zlabel('Z')
ax.set_zlim(-10, 10)

plt.show()