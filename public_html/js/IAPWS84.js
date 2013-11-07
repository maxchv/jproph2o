//# Copyright (c) 2010, Kiran Pashikanti 
//# 
//# Permission to use, copy, modify, and/or distribute this software for any 
//# purpose with or without fee is hereby granted, provided that the above 
//# copyright notice and this permission notice appear in all copies. 
//# 
//# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES 
//# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF 
//# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY 
//# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES 
//# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN 
//# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR 
//# IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 
//# 

//"""
//IAPWS84.py
//
//Implementation of the IAPWS-84 specification for water and steam properties
//"""

var IAPWS84_CRITICAL_TEMPERATURE = 647.096;                          //K
var IAPWS84_CRITICAL_PRESSURE = 22.064;                              //MPa
var IAPWS84_CRITICAL_DENSITY = 322.0;                                //kg/m3
var IAPWS84_ALPHA_0 = 1000.0;                                        //J/kg
var IAPWS84_PHI_0 = IAPWS84_ALPHA_0/IAPWS84_CRITICAL_TEMPERATURE;    //J/kg-K

function vapor_pressure(T) {
    var tau = 1.0 - (T / IAPWS84_CRITICAL_TEMPERATURE);
    var term1 = 0.0;
    term1 += -7.85951783 * tau;
    term1 += 1.84408259 * Math.pow(tau, 1.5);
    term1 += -11.7866497 * Math.pow(tau, 3.0);
    term1 += 22.6807411 * Math.pow(tau, 3.5);
    term1 += -15.9618719 * Math.pow(tau, 4.0);
    term1 += 1.80122502 * Math.pow(tau, 7.5);
    var fact2 = (IAPWS84_CRITICAL_TEMPERATURE / T);
    
    return Math.exp(fact2 * term1) * IAPWS84_CRITICAL_PRESSURE;
}

function d_vapor_pressure_dT(T) {
    var Tc = IAPWS84_CRITICAL_TEMPERATURE;
    var tau = 1.0 - (T / Tc);
    var P = vapor_pressure(T);
    var a1 = -7.85951783;
    var a2 = 1.84408259;
    var a3 = -11.7866497;
    var a4 = 22.6807411;
    var a5 = -15.9618719;
    var a6 = 1.80122502;

    var first = a1 * tau;
    first += a2 * Math.pow(tau, 1.5);
    first += a3 * Math.pow(tau, 3);
    first += a4 * Math.pow(tau, 3.5);
    first += a5 * Math.pow(tau, 4);
    first += a6 * Math.pow(tau, 7.5);

    var first_der = a1;
    first_der += a2 * 1.5 * Math.pow(tau, 0.5);
    first_der += a3 * 3.0 * Math.pow(tau, 2.0);
    first_der += a4 * 3.5 * Math.pow(tau, 2.5);
    first_der += a5 * 4.0 * Math.pow(tau, 3.0);
    first_der += a6 * 7.5 * Math.pow(tau, 6.5);
    first_der *=  (-1.0 / IAPWS84_CRITICAL_TEMPERATURE);

    var second = IAPWS84_CRITICAL_TEMPERATURE / T;
    var second_der = -IAPWS84_CRITICAL_TEMPERATURE / Math.pow(T, 2);

    return P * (first * second_der + second * first_der) * 1e6;
}
    
function saturation_temperature(P) {
    if (P > IAPWS84_CRITICAL_PRESSURE) {
        return null;
    } else if (P === IAPWS84_CRITICAL_PRESSURE) {
        return IAPWS84_CRITICAL_TEMPERATURE;
    }

    var Pr = P / IAPWS84_CRITICAL_PRESSURE;
    //Backward correlation - This is close but not good enough
    var Tr = (0.0383862918967118) * Math.log(Pr) + (-0.251741353110079) * Math.pow(Pr, 2) + (0.435703778366649) * (Pr) + 0.794314236176422;
    var T = Tr * IAPWS84_CRITICAL_TEMPERATURE;
    var converged = False;
    var iterations = 0;
    while (!converged && iterations < 100) {
        if (T > IAPWS84_CRITICAL_TEMPERATURE) {
            T = IAPWS84_CRITICAL_TEMPERATURE - 1.0;
        }
        var f = vapor_pressure(T) - P;
        var T_new = T - (f) / (d_vapor_pressure_dT(T) * 1e-6);
        if (abs(T - T_new) < 1e-10) {
            converged = True;
            break;
        }
        T = T_new;
        iterations += 1;
    }
    if (!converged) {
        return null;
    }

    return T;
}

function density_of_saturated_liquid(T) {
    var tau = 1.0 - (T / IAPWS84_CRITICAL_TEMPERATURE);
    var term1 = 1.0;

    term1 += 1.99274064 * Math.pow(tau, (1.0 / 3.0));
    term1 += 1.09965342 * Math.pow(tau, (2.0 / 3.0));
    term1 += -0.510839303 * Math.pow(tau, (5.0 / 3.0));
    term1 += -1.754934790 * Math.pow(tau, (16.0 / 3.0));
    term1 += -45.5170352 * Math.pow(tau, (43.0 / 3.0));
    term1 += -6.74694450e5 * Math.pow(tau, (110.0 / 3.0));

    return term1 * IAPWS84_CRITICAL_DENSITY;
}

function density_of_saturated_vapor(T) {
    var tau = 1.0 - (T / IAPWS84_CRITICAL_TEMPERATURE);
    var term1 = 0.0;

    term1 += -2.03150240 * Math.pow(tau, (2.0 / 6.0));
    term1 += -2.68302940 * Math.pow(tau, (4.0 / 6.0));
    term1 += -5.38626492 * Math.pow(tau, (8.0 / 6.0));
    term1 += -17.2991605 * Math.pow(tau, (18.0 / 6.0));
    term1 += -44.7586581 * Math.pow(tau, (37.0 / 6.0));
    term1 += -63.9201063 * Math.pow(tau, (71.0 / 6.0));

    return Math.exp(term1) * IAPWS84_CRITICAL_DENSITY;
}

function IAPWS84_Alpha(T) {
    var theta = T / IAPWS84_CRITICAL_TEMPERATURE;
    var term1 = -1135.905627715;

    term1 += -5.65134998e-8 * Math.pow(theta, -19.0);
    term1 += 2690.66631 * theta;
    term1 += 127.287297 * Math.pow(theta, 4.5);
    term1 += -135.003439 * Math.pow(theta, 5.0);
    term1 += 0.981825814 * Math.pow(theta, 54.5);

    return IAPWS84_ALPHA_0 * term1;
}

function IAPWS84_Phi(T) {
    var theta = T / IAPWS84_CRITICAL_TEMPERATURE;
    var term1 = 2319.5246;

    term1 += (19.0 / 20.0) * (-5.65134998e-8) * Math.pow(theta, -20.0);
    term1 += 2690.66631 * Math.log(theta);
    term1 += (9.0 / 7.0) * (127.287297) * Math.pow(theta, 3.5);
    term1 += (5.0 / 4.0) * (-135.003439) * Math.pow(theta, 4.0);
    term1 += (109.0 / 107.0) * (0.981825814) * Math.pow(theta, 53.5);
    
    return IAPWS84_PHI_0 * term1;
}

function specific_enthalpy_of_saturated_liquid(T) {
    var t1 = IAPWS84_Alpha(T);
    t1 += (T/density_of_saturated_liquid(T))*d_vapor_pressure_dT(T);
    return t1;
}

function specific_enthalpy_of_saturated_vapor(T) {
    var t1 = IAPWS84_Alpha(T);
    t1 += (T/density_of_saturated_vapor(T))*d_vapor_pressure_dT(T);
    return t1;
}

function specific_entropy_of_saturated_liquid(T) {
    var t1 = IAPWS84_Phi(T);
    t1 += (1/density_of_saturated_liquid(T))*d_vapor_pressure_dT(T);
    return t1;
}

function specific_entropy_of_saturated_vapor(T) {
    var t1 = IAPWS84_Phi(T);
    t1 += (1/density_of_saturated_vapor(T))*d_vapor_pressure_dT(T);
    return t1;
}

function main() {
    var T = 0.0 + 273.15;
    while (T <= 650.0) {
        //P = vapor_pressure(T);
        var sV = specific_enthalpy_of_saturated_vapor(T)*1e-3;
        var sL = specific_enthalpy_of_saturated_liquid(T)*1e-3;
        console.log(T.toFixed(2), sV.toFixed(5), sL.toFixed(8));
        T += 1;
    }
}

module.exports = {
    specific_entropy_of_saturated_vapor: specific_entropy_of_saturated_vapor
};

//main();
