import numpy as np


def lagrange(x_data, y_data, x):
    result = 0
    for i in range(len(x_data)):
        mult = y_data[i]
        for j in range(len(x_data)):
            if j != i:
                mult *= (x - x_data[j]) / (x_data[i] - x_data[j])
        result += mult
    return result


def newton_coefficients(x_data, y_data):
    coeffs = y_data.copy()
    for i in range(1, len(y_data)):
        for j in range(i, len(y_data)):
            coeffs[j] = (coeffs[j] - coeffs[i - 1]) / (x_data[j] - x_data[i - 1])
    return coeffs


def newton_evaluation(coeffs, x_data, x):
    if len(coeffs) > 1:
        return coeffs[0] + (x - x_data[0]) * newton_evaluation(coeffs[1:], x_data[1:], x)
    return coeffs[0]


def neville(x_data, y_data, x):
    if len(x_data) > 1:
        left = (x - x_data[-1]) * neville(x_data[:-1], y_data[:-1], x)
        right = (x_data[0] - x) * neville(x_data[1:], y_data[1:], x)
        return (left + right) / (x_data[0] - x_data[-1])
    return y_data[0]


xValues = np.array([-1, 2, 4, 6, 8], float)
yValues = np.array([1, 4, 16, 36, 64], float)
# print(lagrange(xValues, yValues, 5))

# print(newton_evaluation(newton_coefficients(xValues, yValues), xValues, 5))

# print(neville(xValues, yValues, 5))

# Neville, switch x and y to get the inverse function of f(x) -> g(x)
# Then, an approximation of the root of f(x) is the result of g(0).
