import Interpolation

def first_center_p(xs, ys, x):
    for i in range(len(xs)):
        if x == xs[i]:
            return (ys[i + 1] - ys[i - 1]) / (2 * (xs[i + 1] - xs[i]))


def first_center_pp(xs, ys, x):
    for i in range(len(xs)):
        if x == xs[i]:
            return (ys[i + 1] - (2 * ys[i]) + ys[i - 1]) / (xs[i + 1] - xs[i]) ** 2


def first_fwd_p(xs, ys, x):
    for i in range(len(xs)):
        if x == xs[i]:
            return (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])


def first_fwd_pp(xs, ys, x):
    for i in range(len(xs)):
        if x == xs[i]:
            return (ys[i + 2] - (2 * ys[i + 1]) + ys[i]) / (xs[i + 1] - xs[i]) ** 2


def second_fwd_p(xs, ys, x):
    for i in range(len(xs)):
        if x == xs[i]:
            return ((-3 * ys[i]) + (4 * ys[i + 1]) - ys[i + 2]) / (2 * (xs[i + 1] - xs[i]))


def second_fwd_pp(xs, ys, x):
    for i in range(len(xs)):
        if x == xs[i]:
            return ((2 * ys[i]) - (5 * ys[i + 1]) + (4 * ys[i + 2]) - ys[i + 3]) / (xs[i + 1] - xs[i]) ** 2


def unequal_step_cubic_spline_first_derivative(xs, ys, x):
    k = Interpolation.cubicspline_curvatures(xs, ys)
    # bracket using Bisection
    i = 0
    ip = len(xs) - 1
    while i + 1 != ip:
        mid = i + ip // 2
        if xs[i] < x < xs[mid]:
            ip = mid
        elif xs[mid] < x < xs[ip]:
            i = mid
    xi = xs[i]
    xip = xs[i + 1]
    yi = ys[i]
    yip = ys[i + 1]
    # Evaluate
    first = (k[i] / 6) * (((3*((x - xip)**2)) / (xi - xip)) - xi + xip)
    second = (k[i + 1] / 6) * (((3*((x - xi)**2)) / (xi - xip)) - xi + xip)
    third = (yi - yip) / (xi - xip)
    return first - second + third
