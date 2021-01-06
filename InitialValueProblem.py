import numpy as np


def euler(Fxy, x, y, xStop, h):
    Xs = [x]
    Ys = [y]
    while abs(x - xStop) > 10**-12:
        h = min(h, xStop - x)
        y += Fxy(x, y) * h
        x += h
        Xs.append(x)
        Ys.append(y)
    return np.array(Xs), np.array(Ys)


"""
def F(x, y):
    return 0.5 * y
x, y = euler(F, 0.0, 1.0, 5.0, 0.01)
def analytical(x):
    return 1.0 * np.exp(0.5 * x)
analytical_y = analytical(x)
e = np.abs(analytical_y - y)
"""


"""
def F_prim(x, y):
    F = np.zeros(2)
    F[0] = y[1]
    F[1] = -4 * y[0]
    return F
results = euler(F_prim, 0.0, np.array([1.0, 0.0]), 5.0, 0.2)
x = results[0]
y = results[1][:, 0]
def analytical(x):
    return 1.0 * np.cos(np.sqrt(4) * x)
analytical_y = analytical(x)
e = np.abs(analytical_y - y)
"""


def RK2(F, x, y, xStop, h):
    X = [x]
    Y = [y]
    hd2 = h / 2
    while abs(xStop - x) > 10 ** -12:
        h = min(h, xStop - x)
        y += h * F(x + hd2, y + (hd2 * F(x, y)))
        x += h
        X.append(x)
        Y.append(y)
    return np.array(X), np.array(Y)


def RK4(F, x, y, xStop, h):
    X = [x]
    Y = [y]
    while abs(xStop - x) > 10 ** -12:
        h = min(h, xStop - x)
        K0 = h * F(x, y)
        K1 = h * F(x + (h / 2.0), y + (K0 / 2.0))
        K2 = h * F(x + (h / 2.0), y + (K1 / 2.0))
        K3 = h * F(x + h, y + K2)
        y += (K0 + (2.0 * K1) + (2.0 * K2) + K3) / 6.0
        x += h
        X.append(x)
        Y.append(y)
    return np.array(X), np.array(Y)



