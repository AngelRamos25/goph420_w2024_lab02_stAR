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

    tol = 0.5e-10
    maxIter = 1000
    eps_r = np.zeros(maxIter)
    eps_r_c = 1
    N = 0

    while eps_r_c > tol and N <= maxIter:
        root = x0 - (f(x0)/dfdx(x0))
        eps_r_c = np.abs((root - x0)/root)
        eps_r[N] = eps_r_c
        x0 = root
        N += 1

    eps_r = eps_r[0:N]
    return (x0, N, eps_r)


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


class ObjectiveFunction:
    def __init__(self,freq) -> None:
        self.freq = freq
    def __call__(self,x):
        rho1 = 1800
        rho2 = 2500
        beta1 = 1900
        beta2 = 3200
        H = 4000
        
        A = (H**2)*( (1/(beta1 ** 2)) -  (1/(beta2 ** 2)) )
        B = rho2/rho1
        U = 2 * (np.pi) * self.freq * x 

        Fx1 = x * np.tan(U)
        Fx2 = B * ( np.sqrt(A - (x ** 2)) )
        Fx = Fx1 - Fx2
        return Fx
    
class ObjectiveFunctionDerivative:
    def __init__(self,freq) -> None:
        self.freq = freq
    def __call__(self,x):
        rho1 = 1800
        rho2 = 2500
        beta1 = 1900
        beta2 = 3200
        H = 4000
        
        A = (H**2)*( (1/(beta1 ** 2)) -  (1/(beta2 ** 2)) )
        B = rho2/rho1
        U = 2 * (np.pi) * self.freq * x 

        Fx1 = np.tan(U) + U * ( (np.tan(U))**2 + 1 )
        Fx2 = (B * x)/(np.sqrt(A - (x**2)))
        Fx = Fx1 + Fx2
        return Fx
