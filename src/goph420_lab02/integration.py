import numpy as np


def root_newton_raphson(
        x0: float,
        f,
        dfdx) -> any:
    """ Newton-Raphson method:
    This function performs a root finding using Newton-Raphson method.

    Inputs:
    ---------------
    x0: Intial guess.
    f: Function f(x).
    dfdx: First derivative of f(x).

    Outputs:
    ---------------
    root: Returns the final estimate of the root.
    N: Returns de number of iterations to convergence.
    eps_r: Returns a one dimensional vector with the approximate relative errors.

    Raises:
    -----------------
    TypeError if f or dfdx are not callable.
    ValueError if x0 is not convertible to float.
    """
    if callable(f) is False:
        raise TypeError(
            "f is not calleable"
        )

    if callable(dfdx) is False:
        raise TypeError(
            "dfdx is not calleable"
        )

    float(x0)

    tol = 0.5e-5
    eps_r = np.zeros((100, 1))
    eps_r_c = 1
    N = -1
    while eps_r_c > tol:
        N += 1
        root = x0 - f(x0)/dfdx(x0)
        eps_r_c = np.abs((root - x0)/root)
        eps_r[N] = eps_r_c
        x0 = root
    print(x0)
    return x0


def root_secant_modified(
        x0: float,
        dx: float,
        f) -> any:
    """ Secant method modified:
    This function performs a root finding using Secant method modified.

    Inputs:
    ---------------
    x0: Intial guess.
    dx: Step size for derivative estimtion.
    f: Function f(x).

    Outputs:
    ---------------
    root: Returns the final estimate of the root.
    N: Returns de number of iterations to convergence.
    eps_r: Returns a one dimensional vector with the approximate relative errors.

    Raises:
    -----------------
    TypeError if f(x) is not callable.
    ValueError if x0 or dx are not convertible to float.
    """
    float(dx)
    float(x0)

    if callable(f) is False:
        raise TypeError(
            "f is not calleable"
        )

    tol = 0.5e-5
    eps_r = np.zeros((100, 1))
    eps_r_c = 1
    N = -1
    while eps_r_c > tol:
        N += 1
        root = x0 - f(x0)*dx/(f(x0 + dx) - f(x0))
        eps_r_c = np.abs((root - x0)/root)
        eps_r[N] = eps_r_c
        x0 = root
    print(x0)
    return x0


class function_Test:

    def __init__(self):
        pass

    def __call__(self, x0):
        f = x0 ** 3 + x0 ** 2 + x0 + 1  # Root = -1
        return f


class function_dfdx_Test:

    def __init__(self):
        pass

    def __call__(self, x0):
        f = 3 * x0 ** 2 + 2 * x0 + 1  # Root = -1
        return f
