import numpy as np
from matplotlib import pyplot as plt
from goph420_lab02 import rootFunctions as rgr


def main():

    beta1 = 1900
    H = 4000
    freq = np.array([0.74,0.75,0.8,0.9,0.95,1,2,3,4,5,6,7,8,9,10]) # Defining frequencys
    k = np.array([0, 1, 2]) # Defining modes
    zetas = np.zeros([len(freq),len(k)])
    cL = np.zeros([len(freq),len(k)])
    lamb = np.zeros([len(freq),len(k)])
    
    # Calculating each root for different frecuencies and modes using Newton Raphson methods:
    a = -1
    
    for i in freq:
        a += 1
        b = -1
        for j in k:
            b += 1
            x0 = (2*j+1)/(4*i) - 1e-4 # Initial guess. 
            FX = rgr.ObjectiveFunction(i)
            FXP = rgr.ObjectiveFunctionDerivative(i)
            root ,N , eps = rgr.root_newton_raphson(x0,FX,FXP)
            zetas[a,b] = root # Saving the roots.
            # Calculating the velocities cL.
            cL[a,b] = np.sqrt( ((beta1 ** 2)*(H ** 2))/ (H**2 - zetas[a,b]*(beta1 ** 2)) )
            lamb[a,b] = cL[a,b]/i
    
    # Plotting f vs zeta for each mode:
    plt.figure(figsize=(20, 6))
    plt.plot(freq,zetas[0:len(freq),0],label = 'mode 0')
    plt.plot(freq,zetas[0:len(freq),1],label = 'mode 1')
    plt.plot(freq,zetas[0:len(freq),2],label = 'mode 2')
    plt.grid(min)
    plt.xlabel('Frequencies [Hz]')
    plt.ylabel(r"$\zeta$ [s]")
    plt.legend()
    plt.savefig(
        'C:/Users/mange/Desktop/UoC/Winter 2024/GOPH_420/goph420_w2024_lab02_stAR/figures/f_vs_z.pdf')
    plt.show()

    # Plotting f vs cL for each mode:
    plt.figure(figsize=(20, 6))
    plt.plot(freq,cL[0:len(freq),0],label = 'mode 0')
    plt.plot(freq,cL[0:len(freq),1],label = 'mode 1')
    plt.plot(freq,cL[0:len(freq),2],label = 'mode 2')
    plt.grid(min)
    plt.xlabel('Frequencies [Hz]')
    plt.ylabel("cL [m/s]")
    plt.legend()
    plt.savefig(
        'C:/Users/mange/Desktop/UoC/Winter 2024/GOPH_420/goph420_w2024_lab02_stAR/figures/f_vs_cL.pdf')
    plt.show()

    # Plotting f vs lambda for each mode:
    plt.figure(figsize=(20, 6))
    plt.plot(freq,lamb[0:len(freq),0],label = 'mode 0')
    plt.plot(freq,lamb[0:len(freq),1],label = 'mode 1')
    plt.plot(freq,lamb[0:len(freq),2],label = 'mode 2')
    plt.grid(min)
    plt.xlabel('Frequencies [Hz]')
    plt.ylabel(r"$\lambda_L$ [m]")
    plt.legend()
    plt.savefig(
        'C:/Users/mange/Desktop/UoC/Winter 2024/GOPH_420/goph420_w2024_lab02_stAR/figures/f_vs_lambda.pdf')
    plt.show()


if __name__ == "__main__":
    main()
