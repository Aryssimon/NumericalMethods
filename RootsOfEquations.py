def rootsearch(f, a, b, dx):
    x = a
    while x < b:
        if (f(x) > 0 > f(x + dx)) or (f(x) < 0 < f(x + dx)):
            return x, x + dx
        x += dx
    return None


def bisection(f, x1, x2, tol=1.0e-9):
    if abs(x2 - x1) <= tol:
        return x1, x2
    x3 = (x1 + x2) / 2
    if (f(x1) < 0 < f(x3)) or (f(x1) > 0 > f(x3)):
        return bisection(f, x1, x3)
    return bisection(f, x3, x2)


def regulafalsi(f, x1, x2, tol=1.0e-9):
    if abs(x2 - x1) <= tol:
        return x1, x2
    x3 = x2 - (f(x2) * ((x2 - x1) / (f(x2) - f(x1))))
    if (f(x1) < 0 < f(x3)) or (f(x1) > 0 > f(x3)):
        return bisection(f, x1, x3)
    return bisection(f, x3, x2)


def secant(f, x1, x2, tol=1.0e-9):
    if abs(x2 - x1) <= tol:
        return x1, x2
    x3 = x2 - (f(x2) * ((x2 - x1) / (f(x2) - f(x1))))
    return secant(f, x2, x3)


def newtonraphson(f, df, a, b, tol=1.0e-9):
    xi = a
    xip = (a + b) / 2
    while abs(xi - xip) > tol:
        xi = xip
        if f(xi) == 0.0:
            return xi
        xip = xip - (f(xip) / df(xip))
    return xi


def newtonraphson_improved(f, df, a, b, tol=1.0e-9):
    xi = a
    xip = (a + b) / 2
    while abs(xi - xip) > tol:
        xi = xip
        if f(xi) == 0.0:
            return xi
        xip = xip - (f(xip) / df(xip))
        if xip < a or xip > b:
            bisection_result = bisection(f, a, b)
            return (bisection_result[0] - bisection_result[1]) / 2
    return xi


def mafunc(x):
   return (3 * x) + 1
def derivativemafunc(x):
   return 3
# print(rootsearch(mafunc, -1, 1, 0.01))
# print(bisection(mafunc, -1, 1))
# print(regulafalsi(mafunc, -1, 1))
# print(secant(mafunc, -1, 1))
# print(newtonraphson(mafunc, derivativemafunc, -1, 1))
# print(newtonraphson_improved(mafunc, derivativemafunc, -1, 1))
