import numpy as np
from matplotlib import pyplot as plt
from goph420_lab02 import integration as Itg

x = np.arange(-100, 100, 1)
f = np.arange(-5, 5, 0.1)
x0 = 2
dx = 0.1
y = Itg.Equation1(-5)
y = y(x)

plt.plot(x, y)
plt.show()

beta1 = 1900

H = 4000

roots = np.zeros([len(f), 1])
cL = np.zeros([len(f), 1])
L = np.zeros([len(f), 1])
a = -1

for i in f:
    print(i)
    a += 1
    gs = Itg.Equation1(i)
    roots[a] = Itg.root_secant_modified(x0, dx, gs)
    # roots[a] = Itg.root_secant_modified(x0, dx, gs)
    cL[a] = 1/np.sqrt(1/(beta1**2) - (roots[a]/H)**2)
    L[a] = cL[a]/i

print(roots)
plt.plot(f, roots)
plt.show()

# if __name__ == "__main__":
#   t, v = np.loadtxt('s_wave_data.txt', float, unpack=True)
#  index = plot_raw_cut(t, v)
# plot_convergence_error(t, v, index)
