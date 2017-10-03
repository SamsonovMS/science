"""
Summation of magnetic fields.
"""

import numpy as np
import matplotlib.pyplot as plt


def quiver_(x, z, U, V, m, n=10, scale=False):
    """
    Helper function.
    :param x:
    :param z:
    :param U:
    :param V:
    :param m:
    :param n:
    :param scale:
    :return:
    """
    if scale:
        plt.quiver(x[::n], z[::n], U[::n, ::n], V[::n, ::n], m[::n, ::n],
                   scale=scale)
    else:
        plt.quiver(x[::n], z[::n], U[::n, ::n], V[::n, ::n], m[::n, ::n],
                   scale=500/n)
    l, r, b, t = plt.axis()
    dx, dy = r - l, t - b
    plt.axis([l - 0.1*dx, r + 0.1*dx, b - 0.1*dy, t + 0.1*dy])
    xgrid, ygrid = np.meshgrid(x, z)
    contour = plt.contour(xgrid, ygrid, m)
    plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
    plt.show()
    #plt.streamplot(x[::n], z[::n], U[::n, ::n], V[::n, ::n])
    #plt.show()


# Shock wave (Harris)
num = [400, 0, 200]
z = np.linspace(-5, 5, num=num[2])
x = np.linspace(-5, 10, num=num[0])
dh = 0.75
bzh = np.tanh(x/dh)
plt.plot(x, bzh)
Uh, Vh = np.meshgrid(x, z)
for i in range(num[2]):
    for j in range(num[0]):
        Uh[i][j] = 0.0
        Vh[i][j] = bzh[j]
mh = np.sqrt(np.power(Uh, 2) + np.power(Vh, 2))
n = 10
scale = 300/n
xgrid, ygrid = np.meshgrid(x, z)

quiver_(x, z, Uh, Vh, mh)

# Plasmoid
d = 0.75
dx1 = 5.0
dx2 = 1.0
bz0 = 1.0 # ?
XO = 7.0 # ?
XM = 12.0 # https://docs.google.com/document/d/1Uebu5muzgRV0MCNWihRWR2x5Dqxgk6k5Sj-HnUkoc_M/edit#bookmark=id.gnczmv5jbq8u
XNL = 15.0 # ?
bz2 = bz0*np.tanh((XM - XO)/dx1)/np.tanh((XNL - XM)/dx2)

bx = np.tanh(z/d)
bz = np.zeros(len(x))
for i, xi in enumerate(x):
    if xi < XM:
        bz[i] = -bz0*np.tanh((xi - XO)/dx1)
    else:
        bz[i] = bz2*np.tanh((xi - XNL)/dx2)

U, V = np.meshgrid(x, z)
for i in range(num[2]):
    for j in range(num[0]):
        U[i][j] = bx[i]
        V[i][j] = bz[j]
m = np.sqrt(np.power(U, 2) + np.power(V, 2))

quiver_(x, z, U, V, m)

# Rotation of plasmoid
a = 0.24
x_axis = -6.0
z_axis = 0.0
Ur, Vr = np.meshgrid(x, z)
for i in range(num[2]): #z
    for j in range(num[0]): #x
        r_point = np.sqrt(np.power(x[j] - x_axis, 2)\
                          + np.power(z[i] - z_axis, 2))
        a_point = np.arctan(float(z[i] - z_axis)/(x[j] - x_axis))
        z_plasm = z_axis + r_point*np.sin(a_point - a)
        b_x = np.tanh(float(z_plasm)/d)
        x_plasm = x_axis + r_point*np.cos(a_point - a)
        if x_plasm < XM:
            b_z = -bz0*np.tanh((x_plasm - XO)/dx1)
        else:
            b_z = bz2*np.tanh((x_plasm - XNL)/dx2)
        r_value = np.sqrt(np.power(b_x, 2) + np.power(b_z, 2))
        a_value = np.arctan(b_z/b_x)
        if b_x < 0:
            a_value = np.pi + a_value
        Ur[i][j] = r_value*np.cos(a_value + a)
        Vr[i][j] = r_value*np.sin(a_value + a)
mr = np.sqrt(np.power(Ur, 2) + np.power(Vr, 2))

quiver_(x, z, Ur, Vr, mr)

# Summation (two cases)
Vs = Vr + Vh
Us = Ur + Uh
ms = np.sqrt(np.power(Us, 2) + np.power(Vs, 2))
quiver_(x, z, Us, Vs, ms, 8)

Vs = -Vr + Vh
Us = -Ur + Uh
ms = np.sqrt(np.power(Us, 2) + np.power(Vs, 2))
quiver_(x, z, Us, Vs, ms, 8, 100)