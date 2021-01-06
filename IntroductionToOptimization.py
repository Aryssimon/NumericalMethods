import numpy as np
import math
import matplotlib.pyplot as plt


def bracket(f, x1, h, c):
    x2 = x1 + h
    fx2 = f(x2)
    if fx2 > f(x1):
        h = -h
        x2 = x1 + h
        fx2 = f(x2)
    for i in range(1000):
        h *= c
        x3 = x2 + h
        fx3 = f(x3)
        if fx3 > fx2:
            return x1, x3
        x1 = x2
        x2 = x3
        fx2 = fx3


"""
def func(x):
    return np.sin(3 * x) + np.cos(x)
print(bracket(func, 1.5, 0.001, 1.1))
x = np.arange(0, 5, 0.01)
y = func(x)
plt.plot(x, y, 'r-')
plt.show()
"""


def golden(f, a, b, tol=1.0e-9):
    nIter = int(np.ceil(-2.078087 * np.log(tol / abs(b - a))))
    R = (-1 + np.sqrt(5)) / 2
    C = 1.0 - R
    x1 = R * a + C * b
    x2 = C * a + R * b
    fx1 = f(x1)
    fx2 = f(x2)
    for i in range(nIter):
        if fx1 > fx2:
            a = x1
            x1 = x2
            fx1 = fx2
            x2 = C * a + R * b
            fx2 = f(x2)
            if abs(x1 - x2) < tol:
                return x2, fx2
        else:
            b = x2
            x2 = x1
            fx2 = fx1
            x1 = R * a + C * b
            fx1 = f(x1)
            if abs(x2 - x1) < tol:
                return x1, fx1
    if fx1 < fx2:
        return x1, fx1
    return x2, fx2
