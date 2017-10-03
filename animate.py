import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import pandas as pd
# my modules
import calc
import convert


n_steps = 500
n_part = 3

a = convert.sol_to_xyz(calc.traces(n_steps=n_steps, n_part=n_part), n_steps, n_part)

t = np.array([np.ones(n_part)*i for i in range(n_steps)]).flatten()
df = pd.DataFrame({'time': t , 'x': a[:, 0], 'y': a[:, 1], 'z': a[:, 2]})


def update_graph(num):
    data = df[df['time']==num]
    graph.set_data(data.x, data.y)
    graph.set_3d_properties(data.z)
    title.set_text('3D Test, step={}'.format(num))
    return title, graph


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('3D Test')

data=df[df['time']==0]
graph, = ax.plot(data.x, data.y, data.z, linestyle="", marker="o")

ani = matplotlib.animation.FuncAnimation(fig, update_graph, n_steps - 1,
                               interval=1, blit=False)

plt.show()