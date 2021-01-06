import numpy as np


def gauss_elimination_multiple_b_vectors(A, b):
    n = len(A)
    m = len(A[0])
    # Forward
    for i in range(n - 1):
        for j in range(i + 1, n):
            multiplier = A[j][i] / A[i][i]
            for k in range(i, m):
                A[j][k] = A[j][k] - (A[i][k] * multiplier)
            for h in range(len(b)):
                b[h][j] = b[h][j] - (b[h][i] * multiplier)
    # Backward
    for i in range(m - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            multiplier = A[j][i] / A[i][i]
            for h in range(len(b)):
                b[h][j] = b[h][j] - (b[h][i] * multiplier)
            A[j][i] = 0
        for h in range(len(b)):
            b[h][i] = b[h][i] / A[i][i]
        A[i][i] = 1

    return A, b


def lu_decomposition(A):
    n = len(A)
    m = len(A[0])
    for i in range(n - 1):
        for j in range(i + 1, n):
            coeff = A[j][i] / A[i][i]
            for k in range(i + 1, m):
                A[j][k] = A[j][k] - (A[i][k] * coeff)
            A[j][i] = coeff
    return A


def lu_solve(LU, b):
    # Forward
    for i in range(1, len(b)):
        sum = 0
        for j in range(i):
            sum += LU[i][j] * b[j]
        b[i] = b[i] - sum
    # Backward
    for i in range(len(b) - 1, -1, -1):
        sum = 0
        for j in range(len(b) - 1, i, -1):
            sum += LU[i][j] * b[j]
        b[i] = (b[i] - sum) / LU[i][i]
    return b


def choleski_decomposition(A):
    n = len(A)
    m = len(A[0])
    for i in range(n):
        for j in range(m):
            if j == 0:
                if i == 0:
                    A[0][0] = np.sqrt(A[0][0])
                else:
                    A[i][0] = A[i][0] / A[0][0]
            elif i == j:
                add = 0
                for k in range(j):
                    add += A[j][k] ** 2
                A[j][j] = np.sqrt(A[j][j] - add)
            else:
                add = 0
                for k in range(j):
                    add += A[i][k] * A[j][k]
                A[i][j] = (A[i][j] - add) / A[j][j]
    return A


def solve_tridiagonal(c, d, e, b):
    # Forward
    for i in range(len(c)):
        b[i + 1] = b[i + 1] - (c[i] * b[i])
    # Backward
    b[-1] = b[-1] / d[-1]
    for i in range(len(b) - 2, -1, -1):
        b[i] = (b[i] - (e[i] * b[i + 1])) / d[i]
    return b


matrix = np.array([[2, 2, 3], [4, 5, 6], [7, 5, 9]], float)
b_matrix = np.array([[15, 32, 44]], float)
# print(gauss_elimination_multiple_b_vectors(matrix, b_matrix))

b_vector = np.array([15, 32, 44], float)
# print(lu_solve(lu_decomposition(matrix), b_vector))

matrix_symmetric_and_positive_definite = np.array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]], float)
# print(choleski_decomposition(matrix_symmetric_and_positive_definite))

c = np.array([-1, 1, -1], float)
d = np.array([1, 2, 1, -1], float)
e = np.array([2, 2, 3], float)
b_vector_tridiag = np.array([8, 12, 22, -9], float)
# print(solve_tridiagonal(c, d, e, b_vector_tridiag))