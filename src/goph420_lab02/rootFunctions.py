import numpy as np


def root_newton_raphson(
        x0: float,
        f,
        dfdx):
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
    eps_r = 0
    eps_r_c = 1

    while eps_r_c > tol:
        root = x0 - f(x0)/dfdx(x0)
        eps_r_c = np.abs((root - x0)/root)
        eps_r = [eps_r, eps_r_c]
        x0 = root
    return x0


def root_secant_modified(
        x0: float,
        dx: float,
        f):
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
    eps_r = 0
    eps_r_c = 1
    while eps_r_c > tol:
        root = x0 - f(x0)*dx/(f(x0 + dx) - f(x0))
        eps_r_c = np.abs((root - x0)/root)
        eps_r = [eps_r, eps_r_c]
        x0 = root
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


class function_Test_2:
    def __init__(self) -> None:
        pass

    def __call__(self, x0) -> float:
        f = x0**3-100*x0**2-x0+100
        return f


class function_dfdx_Test_2:
    def __init__(self) -> None:
        pass

    def __call__(self, x0) -> float:
        f = 3*x0**2 - 200*x0 - 1
        return f


class Equation1:
    def __init__(self, f: float = 0.001):
        self.f = f

    def __call__(self, x) -> float:
        H = 4000
        beta1 = 1900
        beta2 = 3200
        rho1 = 1800
        rho2 = 2500

        tan = np.tan(2*np.pi*self.f*x)
        rhoRatio = (rho2/rho1)
        res = 1/(beta1**2) - 1/(beta2**2)
        gs = x*tan - rhoRatio*np.sqrt((H**2)*res - x**2)
        return gs


class Equation1_derivative:
    def __init__(self, f=0.001) -> None:
        self.f = f

    def __call__(self, x):
        rho1 = 1800
        rho2 = 2500
        beta1 = 1900
        beta2 = 3200
        H = 4000

        rhoRatio = (rho2/rho1)
        res = 1/(beta1**2) - 1/(beta2**2)
        U = 2*np.pi*self.f*x
        tan = np.tan(U)
        gsp = tan + U*(1 + tan ** 2) + rhoRatio * \
            (x / np.sqrt((H**2)*res - x**2))
        return gsp
