import numpy as np


def jacobi(A, b, epsilon):
    X = [np.zeros(len(b), float)]
    for k in range(10000):
        row = np.zeros(len(b), float)
        for i in range(len(row)):
            val = b[i]
            for j in range(len(row)):
                if j != i:
                    val -= A[i][j] * X[k][j]
            row[i] = val / A[i][i]
        X.append(row)
        if np.abs(np.max(X[k + 1] - X[k])) < epsilon:
            return X[k + 1]
    return np.zeros(len(b))


def gauss_seidel_and_sor(A, b, epsilon, w=1):
    X = np.zeros(len(b), float)
    for k in range(10000):
        oldX = X.copy()
        for i in range(len(X)):
            val = b[i]
            for j in range(len(X)):
                if j != i:
                    val -= A[i][j] * X[j]
            X[i] = (1 - w) * X[i] + val * (w / A[i][i])

        if np.abs(np.max(X - oldX)) < epsilon:
            return X
    return np.zeros(len(b))


matrix = np.array([[4, -1, 1], [-1, 4, -2], [1, -2, 4]], float)
b = np.array([12, -1, 5], float)
# print(jacobi(matrix, b, 0.0001))
# print(gauss_seidel_and_sor(matrix, b, 0.0001))