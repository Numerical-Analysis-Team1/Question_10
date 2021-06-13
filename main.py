import random
from sympy import exp
from sympy import log
import sympy as sp
from sympy.utilities.lambdify import lambdify

x = sp.symbols('x')
function = (x * exp(-x) + log(x ** 2)) * ((2 * (x ** 3)) + (2 * (x ** 2)) - 3 * x - 5)


def Newton_Raphson(f, start_point, end_point, epsilon):
    f_prime = f.diff(x)
    f = lambdify(x, f)
    f_prime = lambdify(x, f_prime)
    count = 0
    X = random.uniform(start_point, end_point)
    while abs((X - round((f(X) / f_prime(X)), 4)) - X) > epsilon and f(X) != 0.0 and count < 100:
        X = X - round((f(X) / f_prime(X)), 4)
        print("{0}) {1} ".format(count, X))
        count += 1
    if count >= 100:
        print("Not converge.")
        return False
    elif f(round(X, 3)) != 0:
        print("Number of iteration: ", count)
    return round(X, 4)


def secant_method(f, start_point, end_point, epsilon):
    f = lambdify(x, f)
    count, x_minus1, X = 0, start_point, end_point
    if (f(X) - f(x_minus1)) != 0.0:
        x1 = (x_minus1 * f(X) - X * f(x_minus1)) / (f(X) - f(x_minus1))
    else:
        x1 = 0.0
    while abs(X - x1) > epsilon and count < 100:
        x_minus1 = X
        X = x1
        x1 = (x_minus1 * f(X) - X * f(x_minus1)) / (f(X) - f(x_minus1))
        print("{0}) {1} ".format(count, x1))
        count += 1
    if count >= 100:
        print("Not converge.")
        return False
    elif f(round(x1, 3)) != 0:
        print("Number of iteration : {}".format(count))
    return round(x1, 4)


def ZeroFunction(F, start_point, end_point, epsilon):
    F_prime = F.diff(x)
    f = lambdify(x, F)
    f_prime = lambdify(x, F_prime)

    def byMethod(method):
        numOfZero = 0
        bound = start_point
        while bound <= end_point:
            if f(bound) * f(round(bound + 0.1, 1)) < 0:  # There is a section point
                zero = method(F, bound, round(bound + 0.1, 1), epsilon)
                if zero is not False:
                    print("X" + str(numOfZero) + "=" + formate(zero))
                    numOfZero += 1
                    print()
            if f(bound) * f(round(bound + 0.1, 1)) == 0:
                zero = method(F, bound, round(bound + 0.2, 1), epsilon)
                if zero is not False:
                    print("X" + str(numOfZero) + "=" + formate(zero))
                    numOfZero += 1
                    print()
                bound = round(bound + 0.1, 1)
            bound = round(bound + 0.1, 1)
        # calc derived and check also
        bound = start_point
        while bound <= end_point:
            if f_prime(bound) * f_prime(round(bound + 0.1, 1)) < 0:
                zero = method(F_prime, bound, round(bound + 0.1, 1), epsilon)
                if zero is not False and f(zero) == 0:
                    print("X" + str(numOfZero) + "=" + formate(zero))
                    numOfZero += 1
                    print()
                else:
                    print("No more solution.")
            if f_prime(bound) * f_prime(bound + 0.1) == 0:
                zero = method(F_prime, bound, round(bound + 0.2, 1), epsilon)
                bound = round(bound + 0.1, 1)
                if zero is not False and f(zero) == 0:
                    print("X" + str(numOfZero) + "=" + formate(zero))
                    numOfZero += 1
                    print()
            bound = round(bound + 0.1, 1)

    # NEWTON RAPHSON METHOD
    print("\n*****Newton-Raphson Method*****")
    byMethod(Newton_Raphson)

    # SECANT METHOD
    print("\n*****Secant Method*****")
    byMethod(secant_method)


def Integration_Simpson_Method(f, n, rng):
    f = lambdify(x, f)
    a, b = rng[0], rng[1]
    h = (b - a) / n
    s = f(a) + f(b)
    X = a
    for i in range(0, n - 1):
        X += h
        s += 4 * f(X) if i % 2 == 0 else 2 * f(X)
        print(str((h / 3) * s))
    print()
    print("Integration Value by Simpson Method = " + str(round((h / 3) * s),3))


def formate(num):
    return str(round(num, 3)) + "00000" + "14" + "1101"


ZeroFunction(function, 0.00001, 1.5, 10 ** -4)
Integration_Simpson_Method(function, 10, [0.5, 1])
print(4)
