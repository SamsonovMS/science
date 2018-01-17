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


def shock_wave(x, z, dh):
    """
    Harris 1964 model.
    """
    bz = np.tanh(x / dh)
    u, v = np.meshgrid(x, z)
    for i in range(len(z)):
        for j in range(len(x)):
            u[i][j] = 0.0
            v[i][j] = bz[j]
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    return u, v, m


def current_sheet(x, z, dz):
    """
    Magnetic field of a structure in the solar wind.
    Harris 1964 model.
    """
    bx = np.tanh(z / dz)
    u, v = np.meshgrid(x, z)
    for i in range(len(z)):
        for j in range(len(x)):
            u[i][j] = bx[i]
            v[i][j] = 0.0
    m = np.sqrt(np.power(u, 2) + np.power(v, 2))
    quiver_(x, z, u, v, m, name="ratated plasmoid")
    return u, v, m


def rotate_current_sheet(x, z, dz, a, x_axis, z_axis):
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
            z_cs = z_axis + r_point*np.sin(a_point - a)
            b_x = np.tanh(float(z_cs)/dz)
            b_z = 0.0
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


def main(x, z, case, dh, dz, angle, x_axis, z_axis):
    """
    Calc and show shock wave, plasmoid, rotated plasmoid and its summ with the shock wave in two cases.
    """
    uh, vh, mh = shock_wave(x, z, dh)
    ur, vr, mr = rotate_current_sheet(x, z, dz, angle, x_axis, z_axis)
    u, v, m = summ_fields(uh, vh, ur, vr, case)
    return u, v, m


if __name__ == "__main__":
    print("run fields.py")

    case = 1
    num = [400, 0, 200]
    min_x = -30
    max_x = 60
    min_z = -40
    max_z = 40
    x = np.linspace(min_x, max_x, num=num[0])
    z = np.linspace(min_z, max_z, num=num[2])

    n_steps = 60
    n_part = 5  # > 1
    t_koeff = 0.5  # bigger => more calculations
    delta = 1.0
    theta = 1.0
    dh = 0.75
    dz = 5
    angle = -np.pi/6
    x_axis = -1.0
    z_axis = 0.0

    main(x, z, case, dh, dz, angle, x_axis, z_axis)

