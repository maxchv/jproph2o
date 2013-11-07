#!/usr/bin/env python

# Copyright (c) 2010, Kiran Pashikanti 
# 
# Permission to use, copy, modify, and/or distribute this software for any 
# purpose with or without fee is hereby granted, provided that the above 
# copyright notice and this permission notice appear in all copies. 
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES 
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY 
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES 
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN 
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR 
# IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 
# 

"""
IAPWS84.py

Implementation of the IAPWS-84 specification for water and steam properties
"""

from math import exp, log

IAPWS84_CRITICAL_TEMPERATURE = 647.096                          #K
IAPWS84_CRITICAL_PRESSURE = 22.064                              #MPa
IAPWS84_CRITICAL_DENSITY = 322.0                                #kg/m3
IAPWS84_ALPHA_0 = 1000.0                                        #J/kg
IAPWS84_PHI_0 = IAPWS84_ALPHA_0/IAPWS84_CRITICAL_TEMPERATURE    #J/kg-K

def vapor_pressure(T):
    tau = 1.0 - (T/IAPWS84_CRITICAL_TEMPERATURE)
    term1 = 0.0
    term1 += -7.85951783*tau
    term1 +=  1.84408259*(tau**1.5)
    term1 += -11.7866497*(tau**3.0)
    term1 +=  22.6807411*(tau**3.5)
    term1 += -15.9618719*(tau**4.0)
    term1 +=  1.80122502*(tau**7.5)
    fact2 = (IAPWS84_CRITICAL_TEMPERATURE/T)
    return exp(fact2*term1)*IAPWS84_CRITICAL_PRESSURE

def d_vapor_pressure_dT(T):
    Tc = IAPWS84_CRITICAL_TEMPERATURE
    tau = 1.0 - (T/Tc)    
    P = vapor_pressure(T)
    a1 = -7.85951783
    a2 = 1.84408259
    a3 = -11.7866497
    a4 = 22.6807411
    a5 = -15.9618719
    a6 = 1.80122502

    first = a1*tau
    first += a2*(tau**1.5)
    first += a3*(tau**3)
    first += a4*(tau**3.5)
    first += a5*(tau**4)
    first += a6*(tau**7.5)

    first_der = a1
    first_der += a2*1.5*(tau**0.5)
    first_der += a3*3.0*(tau**2.0)
    first_der += a4*3.5*(tau**2.5)
    first_der += a5*4.0*(tau**3.0)
    first_der += a6*7.5*(tau**6.5)
    first_der = first_der*(-1.0/IAPWS84_CRITICAL_TEMPERATURE)

    second = IAPWS84_CRITICAL_TEMPERATURE/T
    second_der = -IAPWS84_CRITICAL_TEMPERATURE/(T**2)

    return P*(first*second_der + second*first_der)*1e6
    

def saturation_temperature(P):
    if P > IAPWS84_CRITICAL_PRESSURE:
        return None
    elif P == IAPWS84_CRITICAL_PRESSURE:
        return IAPWS84_CRITICAL_TEMPERATURE

    Pr = P/IAPWS84_CRITICAL_PRESSURE
    #Backward correlation - This is close but not good enough
    Tr = (0.0383862918967118)*log(Pr)+(-0.251741353110079)*(Pr**2)+(0.435703778366649)*(Pr)+0.794314236176422
    T = Tr*IAPWS84_CRITICAL_TEMPERATURE
    converged = False
    iterations = 0
    while not converged and iterations < 100:
        if T > IAPWS84_CRITICAL_TEMPERATURE:
            T = IAPWS84_CRITICAL_TEMPERATURE - 1.0
        f = vapor_pressure(T) - P
        T_new = T - (f)/(d_vapor_pressure_dT(T)*1e-6)
        if abs(T-T_new) < 1e-10:
            converged = True
            break
        T = T_new
        iterations += 1

    if converged is False:
        return None
    
    return T


def density_of_saturated_liquid(T):
    tau = 1.0 - (T/IAPWS84_CRITICAL_TEMPERATURE)
    term1 = 1.0
    term1 += 1.99274064*(tau**(1.0/3.0))
    term1 += 1.09965342*(tau**(2.0/3.0))
    term1 += -0.510839303*(tau**(5.0/3.0))
    term1 += -1.754934790*(tau**(16.0/3.0))
    term1 += -45.5170352*(tau**(43.0/3.0))
    term1 += -6.74694450e5*(tau**(110.0/3.0))
    return term1*IAPWS84_CRITICAL_DENSITY

def density_of_saturated_vapor(T):
    tau = 1.0 - (T/IAPWS84_CRITICAL_TEMPERATURE)
    term1 = 0.0
    term1 += -2.03150240*(tau**(2.0/6.0))
    term1 += -2.68302940*(tau**(4.0/6.0))
    term1 += -5.38626492*(tau**(8.0/6.0))
    term1 += -17.2991605*(tau**(18.0/6.0))
    term1 += -44.7586581*(tau**(37.0/6.0))
    term1 += -63.9201063*(tau**(71.0/6.0))
    return exp(term1)*IAPWS84_CRITICAL_DENSITY


def IAPWS84_Alpha(T):
    theta = T/IAPWS84_CRITICAL_TEMPERATURE
    term1 = -1135.905627715
    term1 += -5.65134998e-8*(theta**-19.0)
    term1 += 2690.66631*theta
    term1 += 127.287297*(theta**4.5)
    term1 += -135.003439*(theta**5.0)
    term1 += 0.981825814*(theta**54.5)
    return IAPWS84_ALPHA_0*term1

def IAPWS84_Phi(T):
    theta = T/IAPWS84_CRITICAL_TEMPERATURE
    term1 = 2319.5246
    term1 += (19.0/20.0)*(-5.65134998e-8)*(theta**-20.0)
    term1 += 2690.66631*log(theta)
    term1 += (9.0/7.0)*(127.287297)*(theta**3.5)
    term1 += (5.0/4.0)*(-135.003439)*(theta**4.0)
    term1 += (109.0/107.0)*(0.981825814)*(theta**53.5)
    return IAPWS84_PHI_0*term1

def specific_enthalpy_of_saturated_liquid(T):
    t1 = IAPWS84_Alpha(T) 
    t1 += (T/density_of_saturated_liquid(T))*d_vapor_pressure_dT(T)
    return t1

def specific_enthalpy_of_saturated_vapor(T):
    t1 = IAPWS84_Alpha(T) 
    t1 += (T/density_of_saturated_vapor(T))*d_vapor_pressure_dT(T)
    return t1

def specific_entropy_of_saturated_liquid(T):
    t1 = IAPWS84_Phi(T) 
    t1 += (1/density_of_saturated_liquid(T))*d_vapor_pressure_dT(T)
    return t1

def specific_entropy_of_saturated_vapor(T):
    t1 = IAPWS84_Phi(T) 
    t1 += (1/density_of_saturated_vapor(T))*d_vapor_pressure_dT(T)
    return t1


if __name__ == '__main__':
    T = 0.0 + 273.15
    while T <= 650.0:
        #P = vapor_pressure(T)
        sV = specific_enthalpy_of_saturated_vapor(T)*1e-3
        sL = specific_enthalpy_of_saturated_liquid(T)*1e-3
        print "%.2f %.5f %.8f" % (T, sV, sL)
        T += 1
    
    
    
