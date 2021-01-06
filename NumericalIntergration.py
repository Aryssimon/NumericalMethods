import numpy as np


def trapezoid_composite(f, a, b, n):
    h = (b - a) / n
    add = 0
    fa = f(a)
    for i in range(n):
        fah = f(a + h)
        add += fa + fah
        a += h
        fa = fah
    return add * h / 2


def trapezoid_recursive(f, a, b, k):
    h = (b - a) / (2 ** (k - 1))
    if k == 1:
        return (f(a) + f(b)) * h / 2
    add = 0
    for i in range(1, 2 ** (k - 2) + 1):
        add += f(a + (((2 * i) - 1) * h))
    return (trapezoid_recursive(f, a, b, k - 1) / 2) + (h * add)


def simpson(f, a, b, n):
    # Using Simpspn 1/3
    h = (b - a) / n
    add = f(a) + f(b)
    for i in range(1, n):
        add += (2 + (2 * (i % 2))) * f(a + h)
        a += h
    return (h / 3) * add


def simpson_with_odd_number(f, a, b, n):
    # Integrate the 3 first using Simpson 3/8
    h = (b - a) / n
    threefirst = ((3 * h) / 8) * (f(a) + (3 * f(a + h)) + (3 * f(a + (2 * h))) + f(a + (3 * h)))
    return threefirst + simpson(f, a + (3 * h), b, n - 3)


def romberg(f, a, b, n):
    h = b - a
    R = np.zeros((n, n), float)
    R[0][0] = (h / 2) * (f(a) + f(b))
    for i in range(2, n + 1):
        add = 0
        for j in range(1, 2 ** (i - 2) + 1):
            add += f(a + ((((2 * j) - 1) * h) / (2 ** (i - 1))))
        R[i - 1][0] = (R[i - 2][0] / 2) + (h / (2 ** (i - 1))) * add
    for j in range(1, n):
        for i in range(j, n):
            p = j * 2
            R[i][j] = (((2 ** p) * R[i][j - 1]) - R[i - 1][j - 1]) / ((2 ** p) - 1)
    return R[-1][-1]


def func(x):
    return np.cos(2 / np.cos(x))
# print(trapezoid_composite(func, 0, 1, 10))
# print(trapezoid_recursive(func, 0, 1, 20))
# print(simpson(func, 0, 1, 20))
# print(simpson_with_odd_number(func, 0, 1, 21))
# print(romberg(func, 0, 1, 4))
