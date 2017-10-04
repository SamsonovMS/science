"""
Summation of magnetic fields.
"""

import numpy as np
import matplotlib.pyplot as plt


def quiver_(x, z, u, v, m, n=8, scale=False, name="field"):
    """
    Function for visualisations.
    """
    if scale:
        plt.quiver(x[::n], z[::n], u[::n, ::n], v[::n, ::n], m[::n, ::n],
                   scale=scale)
    else:
        plt.quiver(x[::n], z[::n], u[::n, ::n], v[::n, ::n], m[::n, ::n],
                   scale=500/n)
    l, r, b, t = plt.axis()
    dx, dy = r - l, t - b
    plt.axis([l - 0.1*dx, r + 0.1*dx, b - 0.1*dy, t + 0.1*dy])
    xgrid, ygrid = np.meshgrid(x, z)
    contour = plt.contour(xgrid, ygrid, m)
    plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
    plt.savefig(name + '.png')
    plt.show()
    # plt.streamplot(x[::n], z[::n], U[::n, ::n], V[::n, ::n])
    # plt.show()


def shock_wave(num, x, z, dh):
    """
    Harris 1964 model.
    """
    bzh = np.tanh(x/dh)
    u, v = np.meshgrid(x, z)
    for i in range(num[2]):
        for j in range(num[0]):
            u[i][j] = 0.0
            v[i][j] = bzh[j]
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    return u, v, m


def plasmoid(x, z, d, dx1, dx2, bz0, xo, xm, xnl, bz2):
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
    """
    bx = np.tanh(z/d)
    bz = np.zeros(len(x))
    for i, xi in enumerate(x):
        if xi < xm:
            bz[i] = -bz0*np.tanh((xi - xo)/dx1)
        else:
            bz[i] = bz2*np.tanh((xi - xnl)/dx2)

    u, v = np.meshgrid(x, z)
    for i in range(len(z)):
        for j in range(len(x)):
            u[i][j] = bx[i]
            v[i][j] = bz[j]
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    return u, v, m


def rotate_plasmoid(x, z, d, dx1, dx2, bz0, xo, xm, xnl, bz2, a, x_axis, z_axis):
    """
    Rotation of plasmoid around (x_axis, z_axis) on angle a.
    """
    u, v = np.meshgrid(x, z)
    for i in range(len(z)):
        for j in range(len(x)):
            r_point = np.sqrt(np.power(x[j] - x_axis, 2)\
                              + np.power(z[i] - z_axis, 2))
            a_point = np.arctan(float(z[i] - z_axis)/(x[j] - x_axis))
            if x[j] < x_axis:
                a_point += np.pi
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
    quiver_(x, z, u, v, m, name="ratated plasmoid")
    return u, v, m


def summ_fields(uh, vh, ur, vr, case):
    """
    Two cases.
    """
    if case == 1:
        v = vr + vh
        u = ur + uh
        m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    else:
        v = -vr + vh
        u = -ur + uh
        m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    return u, v, m


def main(num, x, z, case, dh, d, dx1, dx2, bz0, xo, xm, xnl, bz2, angle, x_axis, z_axis):
    """
    Calc and show shock wave, plasmoid, rotated plasmoid and its summ with the shock wave in two cases.
    """
    uh, vh, mh = shock_wave(num, x, z, dh)
    ur, vr, mr = rotate_plasmoid(x, z, d, dx1, dx2, bz0, xo, xm, xnl, bz2, angle, x_axis, z_axis)
    u, v, m = summ_fields(uh, vh, ur, vr, case)
    return u, v, m


if __name__ == "__main__":
    print("run fields.py")
    main()

