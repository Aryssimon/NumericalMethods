import numpy as np
import matplotlib.pyplot as plt
import LinearSystems


def cubicspline_curvatures(x_data, y_data):
    n = len(x_data) - 1
    m = len(x_data)
    c = np.zeros(m - 1, float)
    d = np.ones(m, float)
    e = np.zeros(m - 1, float)
    k = np.ones(m, float)
    for i in range(m - 2):
        c[i] = x_data[i] - x_data[i + 1]
    for i in range(1, m - 1):
        d[i] = 2.0 * (x_data[i - 1] - x_data[i + 1])
        e[i] = x_data[i] - x_data[i + 1]
        left = (y_data[i - 1] - y_data[i]) / (x_data[i - 1] - x_data[i])
        right = (y_data[i] - y_data[i + 1]) / (x_data[i] - x_data[i + 1])
        k[i] = 6.0 * (left - right)
    return LinearSystems.solve_tridiagonal(c, d, e, k)


def cubicspline_evaluation(x_data, y_data, k, x):
    # bracket using Bisection
    a = 0
    b = len(x_data) - 1
    while a + 1 != b:
        mid = a + b // 2
        if x_data[a] < x < x_data[mid]:
            b = mid
        elif x_data[mid] < x < x_data[b]:
            a = mid
    i = a
    xi = x_data[i]
    xip = x_data[i + 1]
    yi = y_data[i]
    yip = y_data[i + 1]
    # Evaluate
    first = (k[i] / 6) * ((((x - xip) ** 3) / (xi - xip)) - ((x - xip) * (xi - xip)))
    second = (k[i + 1] / 6) * ((((x - xi) ** 3) / (xi - xip)) - ((x - xi) * (xi - xip)))
    third = ((yi * (x - xip)) - (yip * (x - xi))) / (xi - xip)
    return first - second + third


xValues = np.array([0, 1, 2], float)
yValues = np.array([0, 2, 1], float)

# curvatures = cubicspline_curvatures(xValues, yValues)
# print(cubicspline_evaluation(xValues, yValues, curvatures, 1.5))


def least_squares_fit(x_data, y_data, degree):
    n = len(x_data)
    A = np.zeros((degree + 1, degree + 1), float)
    b = np.zeros(degree + 1, float)
    # Create matrix and vector
    for i in range((2 * degree) + 1):
        suma = 0
        sumb = 0
        for j in range(n):
            suma += x_data[j] ** i
            sumb += (x_data[j] ** i)*y_data[j]
        for j in range(degree + 1):
            for k in range(degree + 1):
                if j + k == i:
                    if A[j][k] != 0.0:
                        break
                    A[j][k] = suma
                    A[k][j] = suma
        if i <= degree:
            b[i] = sumb

    # Solve matrix
    return np.linalg.solve(A, b)


def plotLeastSF():
    x_data = np.array([1.0, 1.1, 1.8, 2.2, 2.5, 3.5, 3.7, 4.0])
    y_data = np.array([6.008, 5.257, 9.549, 11.098, 15.722, 27.130, 28.828, 33.772])
    x = np.arange(1.0, 4.1, 0.1)
    a_linear = least_squares_fit(x_data, y_data, 1)
    a_quadratic = least_squares_fit(x_data, y_data, 2)
    y_linear = x.copy()
    y_quadratic = x.copy()
    for i in range(len(x)):
        y_linear[i] = a_linear[0] + (a_linear[1] * x[i])
        y_quadratic[i] = a_quadratic[0] + (a_quadratic[1] * x[i]) + (a_quadratic[2] * (x[i]**2))

    # Plot
    plt.plot(x_data, y_data, 'o', x, y_linear, 'r-', x, y_quadratic, 'g-')
    plt.title("Linear and Quadratic LeastSF")
    plt.show()


plotLeastSF()



