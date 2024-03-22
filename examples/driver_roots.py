import numpy as np
from matplotlib import pyplot as plt
from goph420_lab02 import rootFunctions as rgr


def main():
    H = 4000
    beta1 = 1900
    beta2 = 3200
    f = 1
    k = 0  # Mode
    x0 = np.pi*(2*k + 1)/(4*f)
    print(x0)

    f = rgr.Equation1(f)
    fp = rgr.Equation1_derivative(f)
    root, N, eps = rgr.root_newton_raphson(x0, f, fp)


if __name__ == "__main__":
    main()
