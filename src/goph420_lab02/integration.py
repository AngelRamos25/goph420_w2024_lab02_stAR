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
    eps_r = np.zeros(())
    N = 0
    while eps_r > 0.5e-5:
        N += 1
        root = x0 - f(x0)/dfdx(x0)
        eps_r[N] = np.abs((root - x0)/root)
        x0 = root

    return (root, N, eps_r)


class functions_Test:

    def __init__(self):
        pass

    def __call__(self, x0):
        f = x0 ** 3 + x0 ** 2 + x0 + 1  # Root = -1
        return f
