"""
Summation of magnetic fields.
"""

import numpy as np
import matplotlib.pyplot as plt


def quiver_(x, z, U, V, m, n=10, scale=False):
    """
    Helper function for visualisations.
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


def shock_wave(num, x, z, dh=0.75, n=10, scale=300):
    """
    Harris 1964 model.
    :return:
    """
    bzh = np.tanh(x/dh)
    u, v = np.meshgrid(x, z)
    for i in range(num[2]):
        for j in range(num[0]):
            u[i][j] = 0.0
            v[i][j] = bzh[j]
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    return u, v, m


def plasmoid(num, x, z, d, dx1, dx2, bz0, xo, xm, xnl, bz2):
    """
    Magnetic field of a structure in the solar wind.
    :param xm: https://docs.google.com/document/d/1Uebu5muzgRV0MCNWihRWR2x5Dqxgk6k5Sj-HnUkoc_M/edit#bookmark=id.gnczmv5jbq8u
    :param d: 0.75
    :param dx1: 5.0
    :param dx2: 1.0
    :param bz0: 1.0
    :param xo: 7.0
    :param xm: 12.0
    :param xnl: 15.0
    :return:
    """
    bx = np.tanh(z/d)
    bz = np.zeros(len(x))
    for i, xi in enumerate(x):
        if xi < xm:
            bz[i] = -bz0*np.tanh((xi - xo)/dx1)
        else:
            bz[i] = bz2*np.tanh((xi - xnl)/dx2)

    u, v = np.meshgrid(x, z)
    for i in range(num[2]):
        for j in range(num[0]):
            u[i][j] = bx[i]
            v[i][j] = bz[j]
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))

    return u, v, m


def rotate_plasmoid(num, x, z, d, dx1, dx2, bz0, xo, xm, xnl, bz2, a=0.24, x_axis=-6.0, z_axis=0.0):
    """
    Rotation of plasmoid around (x_axis, z_axis) on angle a.
    """
    u, v = np.meshgrid(x, z)
    for i in range(num[2]): #z
        for j in range(num[0]): #x
            r_point = np.sqrt(np.power(x[j] - x_axis, 2)\
                              + np.power(z[i] - z_axis, 2))
            a_point = np.arctan(float(z[i] - z_axis)/(x[j] - x_axis))
            z_plasm = z_axis + r_point*np.sin(a_point - a)
            b_x = np.tanh(float(z_plasm)/d)
            x_plasm = x_axis + r_point*np.cos(a_point - a)
            if x_plasm < xm:
                b_z = -bz0*np.tanh((x_plasm - xo) / dx1)
            else:
                b_z = bz2*np.tanh((x_plasm - xnl) / dx2)
            r_value = np.sqrt(np.power(b_x, 2) + np.power(b_z, 2))
            a_value = np.arctan(b_z/b_x)
            if b_x < 0:
                a_value = np.pi + a_value
            u[i][j] = r_value*np.cos(a_value + a)
            v[i][j] = r_value*np.sin(a_value + a)
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))

    quiver_(x, z, u, v, m)
    return u, v, m


def summ_fields(x, z, uh, vh, ur, vr):
    """
    Two cases.
    """
    v = vr + vh
    u = ur + uh
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    quiver_(x, z, u, v, m, 8)

    v = -vr + vh
    u = -ur + uh
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    quiver_(x, z, u, v, m, 8)
    return u, v, m


def main():
    """
    Calc and show shock wavep, lasmoid, rotated plasmoid and its summ with the shock wave in two cases.
    """
    num = [400, 0, 200]
    x = np.linspace(-5, 10, num=num[0])
    z = np.linspace(-5, 5, num=num[2])

    uh, vh, mh = shock_wave(num, x, z, dh=0.75, n=10, scale=300)

    d = 0.75
    dx1 = 5.0
    dx2 = 1.0
    bz0 = 1.0
    xo = 7.0
    xm = 12.0
    xnl = 15.0
    bz2 = bz0 * np.tanh((xm - xo) / dx1) / np.tanh((xnl - xm) / dx2)

    plasmoid(num, x, z, d, dx1, dx2, bz0, xo, xm, xnl, bz2)
    ur, vr, mr = rotate_plasmoid(num, x, z, d, dx1, dx2, bz0, xo, xm, xnl, bz2, a=0.30, x_axis=-6.0, z_axis=-3.0)
    summ_fields(x, z, uh, vh, ur, vr)


if __name__ == "__main__":
    print("run fields.py")
    main()

