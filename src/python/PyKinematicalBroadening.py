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

#file = 'test_spectrum.spec'
#o = open(file)
#r = o.readlines()
#o.close()

#l = [] ; f = []
#for i in enumerate(r):
#    i_split = i[1].split()
##    if i_split[0] != '#':
#        #print(i_split)
#        l.append(i_split[0]) ; f.append(i_split[1])

#l = np.array(l, dtype=float) ; f = np.array(f, dtype=float) 

# Interpolation for equally spaced wavelength steps
#lambda_o = np.arange(3000,9001,1.)
#fluxes_o = np.interp(lambda_o,l,f)
#Nlambdao = lambda_o.size

#lambda_s = lambda_o
#vc0_gals = 0.0
#Ni_Gauss = 51
#vd_sigma =0.0

#vd = np.arange(0.,1050.,50.)

#n = vd.size
#colors = pl.cm.jet(np.linspace(0,1,n))

#plt.xlim(3000,9000)
#plt.plot(lambda_o,fluxes_o)
#for i in enumerate(vd):
#    Ni_Gauss = max(51,int( i[1] ) + 1)
#    if (Ni_Gauss % 2) == 0:
#        Ni_Gauss += 1
        
#    print('vd = {0:}'.format(i[1]))
#    fluxes_s, IskeepOn = GaussianConvolution( lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma+i[1], Ni_Gauss=Ni_Gauss, fill_val=0.0, verbosity=0 )
#    plt.plot(lambda_s,fluxes_s,label=r'$\sigma_\star$ = {0:}'.format(i), color=colors[i[0]])

#legend = plt.legend(loc=1, prop={'size': 5}, title=r'Velocity dispersion $\sigma_\star$ [km/s]')
#legend.get_title().set_fontsize('5')
#plt.setp( legend.get_title(),fontsize='xx-small')

#plt.xlabel(r'Wavelength [$\AA$]')
#plt.ylabel(r'$F_\lambda$')
