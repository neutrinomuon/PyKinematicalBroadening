#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Revised interface on Sat Jan 30 12:07:21 2021

@author: Jean Gomes

RESUME :  Kinematical broadening of a spectrum in velocity space! Extragalactic Kinematical broadening is a repository for 
applying a kernel in velocity space to models in order to obtain the respective broadened model.

Version: v01

Written: Jean Michel Gomes © Copyright
Created on un Oct 15 10:00:52 WEST 2006
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

import numpy as np
from scipy.special import hermite
from scipy.special import hermitenorm

def broadening_gh(lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma, Ni_Gauss=41, n_hermite=2, coeff_hermite=[0.0, 0.0], fill_val=0.0, verbosity=0):
    
    c  = 2.99792458 * 100000.0
    Pi = np.pi
    
    IsKeepOn = 1

    Nlambdao = np.size(lambda_o)
    Nlambdas = np.size(lambda_s)
    fluxes_s = np.zeros([Nlambdas])

    N_sigmas = 6
    ullambda = -np.float64(N_sigmas)
    uulambda = +np.float64(N_sigmas)
    dulambda = (uulambda - ullambda) / np.float64(Ni_Gauss - 1)

    u_array = np.arange(ullambda,uulambda+dulambda,dulambda)
    coeff_sum = np.array([coeff_hermite[i] * hermitenorm(i)(u_array) for i in range(n_hermite)], dtype=float)
    coeff_sum_total = 0.0
    for i in range(n_hermite):
        coeff_sum_total += coeff_sum[i,:]
    
    # Calculate LOSVD with Gauss-Hermite series
    f_losvd_gauss_hermite = np.exp(-(u_array ** 2 / 2.0)) * (1.0 + coeff_sum_total)
    
    if verbosity == 2:
        print("... DEBUG mode to chek LOSVD shapes")
        f_losvd = np.exp(-(u_array ** 2 / 2.0))
        return u_array,f_losvd,f_losvd_gauss_hermite

    dlambdao = np.gradient(lambda_o)
    v_0_gals = vc0_gals / c

    if v_0_gals != 0.0:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c + v_0_gals)
        lamb_max = lambda_o[Nlambdao - 1] - lambda_o[Nlambdao - 1] * (vd_sigma / c)
    else:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c)
        lamb_max = lambda_o[Nlambdao - 1] - lambda_o[Nlambdao - 1] * (vd_sigma / c + v_0_gals)

    for ilambdas in range(Nlambdas):

        if lamb_min < lambda_s[ilambdas] < lamb_max:

            sumfluxg = 0.0
            # u = ullambda
            
            # Calculate Gauss-Hermite coefficients
            # coeff_sum = sum([coeff_hermite[i] * hermitenorm(i)(u) for i in range(n_hermite)])
            sum_area = 0.0

            # while u <= uulambda:
            for i in range(u_array.size):
                
                # Calculate LOSVD with Gauss-Hermite series
                # f_losvd = np.exp(-(u ** 2 / 2.0)) * (1.0 + coeff_sum_total[])
                
                v = vc0_gals + abs(vd_sigma) * u_array[i]
                l = lambda_s[ilambdas] / (1.0 + v / c)

                if l < lambda_o[0]:
                    f = fluxes_o[0]
                elif l > lambda_o[Nlambdao - 1]:
                    f = fluxes_o[Nlambdao - 1]
                else:
                    N_indice = np.searchsorted(lambda_o, l, side='right') - 1
                    a = (fluxes_o[N_indice + 1] - fluxes_o[N_indice]) / (lambda_o[N_indice + 1] - lambda_o[N_indice])
                    b = (fluxes_o[N_indice] - a * lambda_o[N_indice])
                    f = a * l + b

                sumfluxg = sumfluxg + f * f_losvd_gauss_hermite[i]
                sum_area = sum_area + f_losvd_gauss_hermite[i] * dulambda
                # u = u + dulambda

            fluxes_s[ilambdas] = sumfluxg * dulambda / sum_area #np.sqrt(2.0 * Pi)

        else:
            fluxes_s[ilambdas] = fill_val

    return fluxes_s, IsKeepOn
                
                

def gauss_hermite_convolution(lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma, num_hermite=5, hermite_coeffs=[1.0, 0.0, -0.5, 0.0, 0.375], fill_val=0.0):
    Pi = np.pi
    c = 2.99792458 * 100000.0

    IsKeepOn = 1

    Nlambdao = np.size(lambda_o)
    Nlambdas = np.size(lambda_s)
    fluxes_s = np.zeros([Nlambdas])

    N_sigmas = 6
    ullambda = -np.float64(N_sigmas)
    uulambda = +np.float64(N_sigmas)

    # Compute Gauss-Hermite coefficients
    hermite_coeffs = np.asarray(hermite_coeffs[:num_hermite])
    hermite_coeffs /= np.sqrt(np.sum(hermite_coeffs**2))

    dlambdao = np.gradient(lambda_o)
    v_0_gals = vc0_gals / c

    if v_0_gals != 0.0:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c + v_0_gals)
        lamb_max = lambda_o[Nlambdao - 1] - lambda_o[Nlambdao - 1] * (vd_sigma / c)
    else:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c)
        lamb_max = lambda_o[Nlambdao - 1] - lambda_o[Nlambdao - 1] * (vd_sigma / c + v_0_gals)

    for ilambdas in range(Nlambdas):

        if lamb_min < lambda_s[ilambdas] < lamb_max:

            sumfluxgh = 0.0
            u = ullambda
            while u <= uulambda:

                v = vc0_gals + abs(vd_sigma) * u
                l = lambda_s[ilambdas] / (1.0 + v / c)

                if l < lambda_o[0]:
                    f = fluxes_o[0]
                elif l > lambda_o[Nlambdao - 1]:
                    f = fluxes_o[Nlambdao - 1]
                else:
                    N_indice = np.searchsorted(lambda_o, l, side='right') - 1
                    a = (fluxes_o[N_indice + 1] - fluxes_o[N_indice]) / (lambda_o[N_indice + 1] - lambda_o[N_indice])
                    b = (fluxes_o[N_indice] - a * lambda_o[N_indice])
                    f = a * l + b

                gauss_term = np.exp(-(u ** 2 / 2.0))
                hermite_term = np.sum([hermite(n)(u) * hermite_coeffs[n] for n in range(num_hermite)])
                sumfluxgh += f * gauss_term * hermite_term
                u += 1

            fluxes_s[ilambdas] = sumfluxgh * np.sqrt(2.0) / np.sqrt(Pi)

        else:
            fluxes_s[ilambdas] = fill_val

    return fluxes_s, IsKeepOn

def broadening(lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma, Ni_Gauss=41, fill_val=0.0, verbosity=0):
    Pi = np.pi
    c  = 2.99792458 * 100000.0

    IsKeepOn = 1

    Nlambdao = np.size(lambda_o)
    Nlambdas = np.size(lambda_s)
    fluxes_s = np.zeros([Nlambdas])

    #if Ni_Gauss < vd_sigma / 1.0:
    #    print('[GaussianConvolution] too few Ni_Gauss points => {0:}'.format(Ni_Gauss))
    #    IsKeepOn = 0
    #    return fluxes_s, IsKeepOn

    N_sigmas = 6
    ullambda = -np.float64(N_sigmas)
    uulambda = +np.float64(N_sigmas)
    dulambda = (uulambda - ullambda) / np.float64(Ni_Gauss - 1)

    dlambdao = np.gradient(lambda_o)
    v_0_gals = vc0_gals / c

    if v_0_gals != 0.0:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c + v_0_gals)
        lamb_max = lambda_o[Nlambdao - 1] - lambda_o[Nlambdao - 1] * (vd_sigma / c)
    else:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c)
        lamb_max = lambda_o[Nlambdao - 1] - lambda_o[Nlambdao - 1] * (vd_sigma / c + v_0_gals)

    for ilambdas in range(Nlambdas):

        if lamb_min < lambda_s[ilambdas] < lamb_max:

            sumfluxg = 0.0
            u = ullambda
            while u <= uulambda:

                v = vc0_gals + abs(vd_sigma) * u
                l = lambda_s[ilambdas] / (1.0 + v / c)

                if l < lambda_o[0]:
                    f = fluxes_o[0]
                elif l > lambda_o[Nlambdao - 1]:
                    f = fluxes_o[Nlambdao - 1]
                else:
                    N_indice = np.searchsorted(lambda_o, l, side='right') - 1
                    a = (fluxes_o[N_indice + 1] - fluxes_o[N_indice]) / (lambda_o[N_indice + 1] - lambda_o[N_indice])
                    b = (fluxes_o[N_indice] - a * lambda_o[N_indice])
                    f = a * l + b

                sumfluxg = sumfluxg + f * np.exp(-(u ** 2 / 2.0))
                u = u + dulambda

            fluxes_s[ilambdas] = sumfluxg * dulambda / np.sqrt(2.0 * Pi)

        else:
            fluxes_s[ilambdas] = fill_val

    return fluxes_s, IsKeepOn


def broadening_OLD( lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma, Ni_Gauss=41, fill_val=0.0, verbosity=0 ):
    Pi=3.141592653589793
    c=2.997925*100000.0
    
    IsKeepOn = 1
    
    Nlambdao = np.size(lambda_o)
    Nlambdas = np.size(lambda_s)
    fluxes_s = np.zeros([Nlambdas])
        
    if Ni_Gauss < vd_sigma/1.0:
       print( '[GaussianConvolution] too few Ni_Gauss points => {0:}' .format(Ni_Gauss) )
       IsKeepOn = 0
       return fluxes_s,IsKeepOn
   
    N_sigmas = 6
    ullambda = -np.float64(N_sigmas)
    uulambda = +np.float64(N_sigmas)
    dulambda = (uulambda - ullambda) / np.float64(Ni_Gauss - 1)

    dlambdao = lambda_o[1] - lambda_o[0] # careful ===>> equally spaced 
    v_0_gals = vc0_gals / c
    
    if v_0_gals != 0.0:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c + v_0_gals)
        lamb_max = lambda_o[Nlambdao-1] - lambda_o[Nlambdao-1] * (vd_sigma / c)
    else:
        lamb_min = lambda_o[0] + lambda_o[0] * (vd_sigma / c)
        lamb_max = lambda_o[Nlambdao-1] - lambda_o[Nlambdao-1] * (vd_sigma / c + v_0_gals)

    for ilambdas in range(Nlambdas):

        #print(ilambdas)

        if lambda_s[ilambdas] > lamb_min and lambda_s[ilambdas] < lamb_max:

            sumfluxg = 0.0
            u = ullambda
            while u <= uulambda:

                v  = vc0_gals + abs(vd_sigma) * u
                l  = lambda_s[ilambdas] / (1.0 + v / c)

                if l < lambda_o[0]:
                    f = fluxes_o[0]
                elif l > lambda_o[Nlambdao-1]:
                    f = fluxes_o[Nlambdao-1]
                else:
                    N_indice = int( (l - lambda_o[0] / dlambdao ) )  + 1
                    #print(N_indice)
                    
                    a = (fluxes_o[N_indice] - fluxes_o[N_indice-1]) / (lambda_o[N_indice] - lambda_o[N_indice-1])
                    b  = (fluxes_o[N_indice-1] - a * lambda_o[N_indice-1])
                    f  = a * l + b
                    
                sumfluxg = sumfluxg + f * np.exp(-(u**2 / 2.0))
                u = u + dulambda
            
            fluxes_s[ilambdas] = sumfluxg * dulambda / np.sqrt(2.0 * Pi)

    else:
           fluxes_s[ilambdas] = fill_val

    return  fluxes_s,IsKeepOn

# file = 'test_spectrum.spec'
# o = open(file)
# r = o.readlines()
# o.close()

# l = [] ; f = []
# for i in enumerate(r):
#     i_split = i[1].split()
#     if i_split[0] != '#':
#         # print(i_split)
#         l.append(i_split[0]) ; f.append(i_split[1])

# l = np.array(l, dtype=float) ; f = np.array(f, dtype=float) 

# # Interpolation for equally spaced wavelength steps
# lambda_o = np.arange(3000,9001,1.)
# fluxes_o = np.interp(lambda_o,l,f)
# Nlambdao = lambda_o.size

# lambda_s = lambda_o
# vc0_gals = 0.0
# Ni_Gauss = 11
# vd_sigma = 0.0

# vd = np.arange(0.,550.,50.)

# n = vd.size
# colors = pl.cm.jet(np.linspace(0,1,n))

# plt.xlim(3000,9000)
# plt.plot(lambda_o,fluxes_o)
# for i in enumerate(vd):
#     # Ni_Gauss = max(51,int( i[1] ) + 1)
#     if (Ni_Gauss % 2) == 0:
#         Ni_Gauss += 1
        
#     print('vd = {0:}'.format(i[1]))
#     # fluxes_s, IskeepOn = broadening( lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma+i[1], Ni_Gauss=Ni_Gauss, fill_val=0.0, verbosity=0 )
#     # fluxes_s, IskeepOn = gauss_hermite_convolution( lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma+i[1], num_hermite=5, hermite_coeffs=[1.0, 0.0, -0.5, 0.0, 0.375], fill_val=0.0)
    
    
#     fluxes_s, IskeepOn = broadening_gh(lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma+i[1], Ni_Gauss=Ni_Gauss, n_hermite=3, coeff_hermite=[0.04, 0.04, 0.04], fill_val=0.0, verbosity=0)
        
    
#     plt.plot(lambda_s,fluxes_s,label=r'$\sigma_\star$ = {0:}'.format(i), color=colors[i[0]])

# legend = plt.legend(loc=1, prop={'size': 5}, title=r'Velocity dispersion $\sigma_\star$ [km/s]')
# legend.get_title().set_fontsize('5')
# plt.setp( legend.get_title(),fontsize='xx-small')

# plt.xlabel(r'Wavelength [$\AA$]')
# plt.ylabel(r'$F_\lambda$')
