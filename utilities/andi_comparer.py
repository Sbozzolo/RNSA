#!/usr/bin/env python3


# Compare universal relation with Bauswein, Stergioulas 2017

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def comparer(eos, mass, baus_thres, sim_thres = False, radius = 13.169):

    # Prameters
    R_ref = 13.169

    a = 4.041
    b = 4.568 + (radius - R_ref)*0.1
    c = 0.301
    d = 0.127
    e = 0.026
    M_max = mass
    baus_k = baus_thres/M_max

    # The two kappa functions

    def kappa(jmoment):
        return b/a*1/M_max + jmoment*(M_max)/a

    def kappa2(jmoment):
        return 1 + c*(jmoment**2) - d*(jmoment**4) + e*(jmoment**6)

    # The equation

    func = lambda jmoment : 1 + c*(jmoment)**2 - d*(jmoment)**4 + e*(jmoment)**6 - b/a* 1/M_max - jmoment/a * M_max

    # Plotting

    # xp = np.linspace(0, 2, 100)
    # yp1 = kappa(xp)
    # yp2 = kappa2(xp)

    # p1, = plt.plot(xp, yp1, "orange", label = "Bauswein")
    # p2, = plt.plot(xp, yp2, "teal", label = "Bozzola")
    # plt.xlabel("$\\frac{cJ}{G{M^\\star}^2}$")
    # plt.ylabel("$\\frac{M}{M^\\star}$")
    # plt.title("Threshold Mass for EOS" + eos)
    # plt.legend(handles=((p1,p2)))
    # plt.show()

    # Solve numerically the equation, first try

    initial_guess = 0.5
    jmoment = fsolve(func, initial_guess)

    baus_k = baus_thres/M_max
    sim_k  = sim_thres/M_max

    error1 = np.abs(kappa(jmoment[0]) - baus_k)/np.average([kappa(jmoment[0]), baus_k])

    # Solve numerically the equation, second try

    # initial_guess = 2
    # jmoment = fsolve(func, initial_guess)

    # error2 = np.abs(kappa(jmoment[0]) - baus_k)/np.average([kappa(jmoment[0]), baus_k])

    # if (error1 < error2):
    #     initial_guess = 0.5
    # else:
    #     initial_guess = 2

    jmoment = fsolve(func, initial_guess)

    print("EOS", eos)
    print("Bozzola's value of k: {:.3f}".format(kappa(jmoment[0])))
    print("Bauswein's value of k: {:.3f}".format(baus_k))
    print("Relative error: {:1.2f} %".format(np.abs(kappa(jmoment[0]) - baus_k)/np.average([kappa(jmoment[0]), baus_k])*100))
    if (not sim_thres == False):
        print("Simulations's value of k: {:.3f}".format(sim_k))
        print("Relative error: {:1.2f} %".format(np.abs(kappa(jmoment[0]) - sim_k)/np.average([kappa(jmoment[0]), sim_k])*100))

    return kappa(jmoment[0])

def comp(M_max, R_max):
    c     = 2.998e+10
    G     = 6.674e-8
    M_sun = 1.988e+33
    return G*M_sun/(c*c*1e5)*M_max/R_max

if __name__ == "__main__":
    kappas = [
        comparer("DD2", 2.42, 3.24, 3.35, 13.169),
        comparer("LS220", 2.04, 2.94, 3.05, 12.471),
        comparer("LS375", 2.71, 3.39, 3.65, 13.631),
        comparer("NL3", 2.79, 3.58, 3.85, 14.708),
        comparer("SFHO", 2.06, 2.86, 2.95, 11.750),
        comparer("SFHX", 2.13, 2.95, 3.05, 12.001),
        comparer("TM1", 2.21, 3.25, 3.45, 14.320),
        comparer("TMA", 2.02, 3.08),
        comparer("APR", 2.19, 2.77),
        comparer("SLy4", 2.05, 2.81),
        comparer("ppAPR3", 2.38, 3.00),
        comparer("ppENG", 2.25, 2.93),
        comparer("ppH4", 2.02, 3.08),
        comparer("ppMPA1", 2.47, 3.15),
        comparer("ppMS1", 2.77, 3.59),
        comparer("ppMS1b", 2.76, 3.55),
        comparer("ppEoSa", 2.05, 3.17),
        comparer("ppEosB", 2.35, 3.32),
    ]
    kappas_baus = [
        1.339,
        1.441,
        1.251,
        1.283,
        1.388,
        1.385,
        1.471,
        1.525,
        1.265,
        1.371,
        1.261,
        1.302,
        1.525,
        1.275,
        1.296,
        1.286,
        1.546,
        1.413
    ]
    comps = [
        comp(2.42, 11.87),
        comp(2.04, 10.61),
        comp(2.71, 12.30),
        comp(2.79, 13.39),
        comp(2.06, 10.30),
        comp(2.13, 10.78),
        comp(2.21, 12.50),
        comp(2.02, 12.11),
        comp(2.19, 9.90),
        comp(2.05, 9.97),
        comp(2.38, 10.73),
        comp(2.25, 10.40),
        comp(2.02, 11.72),
        comp(2.47, 11.34),
        comp(2.77, 13.37),
        comp(2.76, 13.28),
        comp(2.05, 12.43),
        comp(2.35, 12.57),
    ]
    p = np.polyfit(comps, kappas, 1)
    print(p)
    print(comps)

    xx = np.linspace(np.min(comps), np.max(comps), 100)
    yy = [-3.49*x + 2.40 for x in xx]

    xx2 = np.linspace(np.min(comps), np.max(comps), 100)
    yy2 = [-3.95*x + 2.50 for x in xx2]

    p1, = plt.plot(comps, kappas, 'o', color = 'red', label = 'Bozzola')
    p2, = plt.plot(comps, kappas_baus, 'o', color = 'blue', label = 'Bauswein')
    plt.plot(xx, yy, 'red')
    plt.plot(xx2, yy2, 'blue')
    plt.legend(handles=((p1,p2)))
    plt.show()
