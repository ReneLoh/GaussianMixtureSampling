import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pylab import figure, cm
from matplotlib.ticker import LogLocator
import csv
from cycler import cycler
import time

sqrt_PI_inv = 1/np.sqrt(2*np.pi)

def gaussmix(mu1,mu2,x,sigs,a_s):
    g1 = sqrt_PI_inv * a_s[0]/sigs[0] * np.exp(-(x-mu1)**2/(2*sigs[0]**2))
    g2 = sqrt_PI_inv * a_s[1]/sigs[1] * np.exp(-(x-mu2)**2/(2*sigs[1]**2))
    return g1+g2

def prior_gauss(mu1,mu2,mu1_mean=0,mu2_mean=0,sigma0=5):
    return 1/(2*np.pi*sigma0**2) * np.exp( -1*( (mu1-mu1_mean)**2+(mu2-mu2_mean)**2 )/(2*sigma0**2) )

def posterior_plot(mu1,mu2, X, likelihood_func, likelihood_params, prior_func, prior_params):
    
    P_prior = prior_func(mu1,mu2, *prior_params)
    # Ls = np.ones((mu1.shape[0], mu2.shape[1]), float)
    Ls = 0
    for x in X:
        Ls += np.log(likelihood_func(mu1,mu2, x, *likelihood_params))
        
    return np.log(P_prior) + Ls

def posterior(mu1,mu2, X, likelihood_func, likelihood_params, prior_func, prior_params):
    
    P_prior = prior_func(mu1,mu2, *prior_params)
    # Ls = np.ones((mu1.shape[0], mu2.shape[1]), float)
    Ls = np.sum(np.log(gaussmix(mu1,mu2, X, *likelihood_params)))
        
    return np.log(P_prior) + Ls


def noisy_force(theta, T, likelihood_func, sigs, a_s, X):
    exp_term1 = np.exp( -(X-theta[0])**2 / (2*sigs[0]**2) ) * (X-theta[0])/sigs[0]**2
    exp_term2 = np.exp( -(X-theta[1])**2 / (2*sigs[1]**2) ) * (X-theta[1])/sigs[1]**2
    P_inv = 1/likelihood_func(*theta, X, sigs, a_s)
    
    sum1 = np.sum(P_inv * exp_term1)
    sum2 = np.sum(P_inv * exp_term2)
        
    return T*sqrt_PI_inv * np.array([ a_s[0]/sigs[0] * sum1, a_s[1]/sigs[1] * sum2 ])


def create_data(mus, sigs, a_s, N):
    # Draw N 1d data points from passed density function of given arguments.
    
    acc_as = np.array( [ sum(a_s[:i]) for i in range(1, len(mus)+1) ] )
    X = []
    for i in range(0,N):
        # choose Gaussian to sample from
        r = np.random.uniform()
        for j in range(0,len(acc_as)):
            if r <= acc_as[j]:
                k = j
                break
        # sample from kth Gaussian
        X += [np.random.normal(mus[k], sigs[k])]
    
    return np.array(X)

def write_to_file(Xdata, filename):
    with open(filename, mode='w') as file:
        writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for x in Xdata:
            writer.writerow([x])
    file.close()




# def OBABO_simu(theta0, p0, N, h, T, gamma, force_func, force_params):
#     start_time = time.time()

#     theta_traj = [theta0]
#     a = np.exp(-gamma*h)
#     theta = np.copy(theta0)
#     print("Theta: ", theta)
#     p = np.copy(p0)
#     print("p: ", p)
#     force = force_func(theta, *force_params)
#     print("Force: ", force)
#     for i in range(0, N):

#         Rn = np.random.normal(size=p.shape)
#         print("Rn: ", Rn)
#         p = np.sqrt(a)*p + np.sqrt((1-a)*T)*Rn + h/2 * force
#         print("p: ", p)
#         theta += h*p
#         print("Theta: ", theta)
#         force = force_func(theta, *force_params)
#         print("Force: ", force)
#         Rn = np.random.normal(size=p.shape)
#         print("Rn: ", Rn)
#         p = p + h/2 * force
#         print("p: ", p)
#         p = np.sqrt(a)*p + np.sqrt((1-a)*T)*Rn
#         print("p: ", p)
        
#         theta_traj += [np.copy(theta)]
#         print("Traj: ", theta_traj)

#         if i % 10000 == 0:
#             print("iteration {} done!\n".format(i))

#     print("Execution took {} minutes!".format((time.time()-start_time)/60))
#     return np.array(theta_traj)

def OBABO_simu(theta0, p0, N, h, T, gamma, force_func, force_params):
    start_time = time.time()

    a = np.exp(-gamma*h)
    theta = np.copy(theta0)
    p = np.copy(p0)
    force = force_func(theta, *force_params)
    
    avg_sum = np.copy(theta)
    T_AVG = [np.copy(avg_sum)]
    for i in range(0, N):

        Rn = np.random.normal(size=p.shape)
        p = np.sqrt(a)*p + np.sqrt((1-a)*T)*Rn + h/2 * force
        theta += h*p
        force = force_func(theta, *force_params)
        Rn = np.random.normal(size=p.shape)
        p = p + h/2 * force
        p = np.sqrt(a)*p + np.sqrt((1-a)*T)*Rn
        
        avg_sum += np.sort(theta)
        T_AVG += [np.copy(avg_sum)]

        if i % 10000 == 0:
            print("iteration {} done!\n".format(i))

    print("Execution took {} minutes!".format((time.time()-start_time)/60))
    T_AVG = [T_AVG[i]/(i+1) for i in range(0,len(T_AVG))]
    return np.array(T_AVG)


def OBABO_simu_MH(theta0, p0, N, h, T, gamma, L, SF, force_func, force_params, 
                  posterior, Xdata, likelihood_func, likelihood_params):
    start_time = time.time()

    a = np.exp(-gamma*h)
    theta_curr = np.copy(theta0)
    p_curr = np.copy(p0)
    force = force_func(theta_curr, *force_params)

    avg_sum = np.copy(theta_curr)
    T_AVG = [np.copy(avg_sum)]
    ctr = 0  # counts acceptances
    
    for i in range(0, N):

        theta = np.copy(theta_curr)
        p = np.copy(p_curr)
        # stores kinetic energies as used by the MH multistep criterion (see Alonso)
        kin_energy = 0
        for k in range(0, L):
            # Perform L OBABO steps
            Rn = np.random.normal(size=p.shape)
            p = np.sqrt(a)*p + np.sqrt((1-a)*T)*Rn  # O
            kin_energy -= 0.5*np.sum(p**2)
            p += h/2 * force   # B
            theta += h*p   # A
            force = force_func(theta, *force_params)
            Rn = np.random.normal(size=p.shape)
            p += h/2 * force  # B
            kin_energy += 0.5*np.sum(p**2)
            p = np.sqrt(a)*p + np.sqrt((1-a)*T)*Rn  # O

        if theta[0]>theta[1]:
            theta = np.flipud(theta)
            p = np.flipud(p)
            
        # MH criterion
        # MH = np.exp(-(1/T) * (U_pot(theta) - U_pot(theta_curr) + kin_energy))
        MH = np.exp(-(1/T) * (-T*posterior(*theta, Xdata, likelihood_func, likelihood_params) 
                              + T*posterior(*theta_curr, Xdata, likelihood_func, likelihood_params) 
                              + kin_energy))        ## is this correct?

        if np.random.rand() < min([1, MH]):
            avg_sum += np.copy(theta)  # accept proposal
            T_AVG += [np.copy(avg_sum)]
            theta_curr = np.copy(theta)  
            if SF==True: 
                p_curr = -1*np.copy(p)
            else:
                p_curr = np.copy(p)
                
            force = force_func(theta_curr, *force_params)
            ctr += 1
        else:
            force = force_func(theta_curr, *force_params)
            avg_sum += np.sort(theta_curr)  
            T_AVG += [np.copy(avg_sum)]  # reject

        if i % 10000 == 0:
            print("iteration {} done!\n".format(i))

    print("acceptance probability was {}".format(ctr/N))
    print("Execution took {} minutes!".format((time.time()-start_time)/60))
    T_AVG = [T_AVG[i]/(i+1) for i in range(0,len(T_AVG))]
    return np.array(T_AVG)
