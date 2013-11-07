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


from google.appengine.ext import webapp
from google.appengine.ext.webapp import util
import logging
from django.utils import simplejson
import IAPWS95
import traceback
import logging

class MainHandler(webapp.RequestHandler):
    def get(self):
        f = file('proph2o.html')
        self.response.out.write(f.read())


class CalcHandler(webapp.RequestHandler):
    def post(self):
        spec1_name = str(self.request.get('specification1_name'))
        spec1_value = float(self.request.get('specification1_value'))
        spec1_unit = str(self.request.get('specification1_unit'))

        spec2_name = str(self.request.get('specification2_name'))
        spec2_value = float(self.request.get('specification2_value'))
        spec2_unit = str(self.request.get('specification2_unit'))
        
        results = function_lookup(spec1_name, spec1_value, spec1_unit,
                                  spec2_name, spec2_value, spec2_unit)
        
        self.response.out.write(simplejson.dumps(results))


def function_lookup(spec1_name, spec1_value, spec1_unit,
                    spec2_name, spec2_value, spec2_unit):

    t = [spec1_name, spec2_name]
    t.sort()
    t = tuple(t)

    MAP = {('Density', 'Temperature'): calculate_Density_Temperature,
           ('Quality', 'Temperature'): calculate_Quality_Temperature,
           ('Pressure', 'Quality'): calculate_Pressure_Quality,
           ('Density', 'Pressure'): calculate_Density_Pressure,
           ('Pressure', 'Temperature'): calculate_Pressure_Temperature,
           
           ('Enthalpy', 'Pressure'): calculate_Enthalpy_Pressure,
           ('Enthalpy', 'Temperature'): calculate_Enthalpy_Temperature,
           ('Enthalpy', 'Quality'): calculate_Enthalpy_Quality,

           ('Entropy', 'Pressure'): calculate_Entropy_Pressure,
           ('Entropy', 'Temperature'): calculate_Entropy_Temperature,
           ('Entropy', 'Quality'): calculate_Entropy_Quality,
           }

    f = MAP.get(t, None)
    if f is None:
        r = new_results()
        r['status'] = 'Those specifications are not supported yet'
    else:
        kw = {spec1_name: (spec1_value, spec1_unit), spec2_name: (spec2_value, spec2_unit)} 
        try:
            r = f(kw)
            
        except Exception, e:
            outs = traceback.format_exc()
            r = new_results()
            r['status'] = 'An internal error has occured, please report this error' + '<pre>' + outs + '</pre>'
            logging.error(r['status'] + str(kw))
            
    return r


def new_results():
    return {"status": None,
            "Quality": ["---", "---", "---"],
            "Temperature": ["---", "---", "---"],
            "Pressure": ["---", "---", "---"],
            "Density": ["---", "---", "---"],
            "Internal_Energy": ["---", "---", "---"],
            "Entropy" : ["---", "---", "---"],
            "Enthalpy": ["---", "---", "---"],
            "Isochoric_heat_capacity": ["---", "---", "---"],
            "Isobaric_heat_capacity": ["---", "---", "---"],
            "Speed_of_sound": ["---", "---", "---"],
            "Joule_Thomson_coefficient": ["---", "---", "---"],
            "Isothermal_throttling_coefficient": ["---", "---", "---"],
            "Isentropic_temperature_pressure_coefficient": ["---", "---", "---"],
            "Second_virial_coefficient": ["---", "---", "---"],
            "Third_virial_coefficient": ["---", "---", "---"]}


def calculate_phase(phase, kw, r):
    rho_value, rho_unit = kw.get('Density')
    T_value, T_unit = kw.get('Temperature')
    try:
        P_value = IAPWS95.pressure(rho_value, T_value)
        U_value = IAPWS95.internal_energy(rho_value, T_value)
        S_value = IAPWS95.entropy(rho_value, T_value)
        H_value = IAPWS95.enthalpy(rho_value, T_value)
        Cv_value = IAPWS95.isochoric_heat_capacity(rho_value, T_value)
        Cp_value = IAPWS95.isobaric_heat_capacity(rho_value, T_value)
        w_value = IAPWS95.speed_of_sound(rho_value, T_value)
        jt_value = IAPWS95.joule_thompson_coefficient(rho_value, T_value)
        it_value = IAPWS95.isothermal_throttling_coefficient(rho_value, T_value)
        iStp_value = IAPWS95.isentropic_temperature_pressure_coefficient(rho_value, T_value)
    except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError):
        r["status"] = 'Unfeasible conditions likely specified, please check input'
        return
    
    r["Temperature"][phase] = T_value
    r["Pressure"][phase] = P_value
    r["Density"][phase] = rho_value
    r["Internal_Energy"][phase] = U_value
    r["Entropy"][phase] = S_value
    r["Enthalpy"][phase] = H_value
    r["Isochoric_heat_capacity"][phase] = Cv_value
    r["Isobaric_heat_capacity"][phase] = Cp_value
    r["Speed_of_sound"][phase] = w_value
    r["Joule_Thomson_coefficient"][phase] = jt_value
    r["Isothermal_throttling_coefficient"][phase] = it_value
    r["Isentropic_temperature_pressure_coefficient"][phase] = iStp_value
    r["Second_virial_coefficient"] = ["---", "---", "---"]
    r["Third_virial_coefficient"] = ["---", "---", "---"]

    return r

def calculate_overall_phase(r):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2

    quality = r["Quality"][VAPOR]
    r["Temperature"][OVERALL] = r["Temperature"][VAPOR]
    r["Pressure"][OVERALL] = r["Pressure"][VAPOR]
    r["Internal_Energy"][OVERALL] = quality*r["Internal_Energy"][VAPOR] + (1.0 - quality)*r["Internal_Energy"][LIQUID]
    r["Entropy"][OVERALL] = quality*r["Entropy"][VAPOR] + (1.0 - quality)*r["Entropy"][LIQUID]
    r["Enthalpy"][OVERALL] = quality*r["Enthalpy"][VAPOR] + (1.0 - quality)*r["Enthalpy"][LIQUID]
    r["Quality"][OVERALL] = 1.0
    
    return r

def calculate_Density_Temperature(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2

    r = new_results()
    T_value, T_unit = kw['Temperature']
    rho_value, rho_unit = kw['Density']

    if str(T_unit) != 'K' or str(rho_unit) != 'kg/m3':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if T_value <= 0.0 or rho_value <= 0.0:
        r['status'] = 'Both temperature and density must have positive values.'
        return r
        
    P_value = IAPWS95.pressure(rho_value, T_value)
    results = IAPWS95.saturation_pressure_at_fixed_temperature(T_value)
    if results is None:
        return calculate_phase(OVERALL, kw, r)

    P_value_sat = results[1]
    rhoV_value = results[2]
    rhoL_value = results[3]
    
    if P_value < P_value_sat:
        r['Quality'][VAPOR] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(VAPOR, kw, r)
    elif P_value > P_value_sat:
        r['Quality'][LIQUID] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(LIQUID, kw, r)
    else:
        kw['Density'] = (rhoV_value, 'kg/m3') 
        r = calculate_phase(VAPOR, kw, r)

        kw['Density'] = (rhoL_value, 'kg/m3')
        r = calculate_phase(LIQUID, kw, r)
        r['status'] = 'Results available'        
        return r


def calculate_Quality_Temperature(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2

    r = new_results()
    T_value, T_unit = kw['Temperature']
    x_value, x_unit = kw['Quality']
    
    if str(T_unit) != 'K' or str(x_unit) != 'Unitless':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if x_value < 0.0 or x_value > 1.0:
        r['status'] = 'Quality must be between 0 and 1'
        return r

    if T_value < 0.0:
        r['status'] = 'Both temperature and density must have positive values.'
        return r

    results = IAPWS95.temperature_quality(T_value, x_value)

    if results is None:
        r['status'] = 'Cannot find solution with those specifications.'
        return r
    
    T_value, P_value, rhoV_value, rhoL_value = results

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoV_value, 'kg/m3')
    r['Quality'][VAPOR] = x_value
    r = calculate_phase(VAPOR, kw, r)

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoL_value, 'kg/m3')
    r['Quality'][LIQUID] = 1.0 - x_value
    r = calculate_phase(LIQUID, kw, r)
    r = calculate_overall_phase(r)

    return r


def calculate_Pressure_Quality(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2

    r = new_results()
    P_value, P_unit = kw['Pressure']
    x_value, x_unit = kw['Quality']
    
    if str(P_unit) != 'kPa' or str(x_unit) != 'Unitless':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if x_value < 0.0 or x_value > 1.0:
        r['status'] = 'Quality must be between 0 and 1'
        return r

    results = IAPWS95.pressure_quality(P_value, x_value)

    if results is None:
        r['status'] = 'Cannot find solution with those specifications.'
        return r
    
    T_value, P_value, rhoV_value, rhoL_value = results

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoV_value, 'kg/m3')
    r['Quality'][VAPOR] = x_value
    r = calculate_phase(VAPOR, kw, r)

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoL_value, 'kg/m3')
    r['Quality'][LIQUID] = 1.0 - x_value
    r = calculate_phase(LIQUID, kw, r)
    r = calculate_overall_phase(r)
    
    return r    


def calculate_Density_Pressure(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2

    r = new_results()
    P_value, P_unit = kw['Pressure']
    rho_value, rho_unit = kw['Density']

    if str(P_unit) != 'kPa' or str(rho_unit) != 'kg/m3':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if P_value <= 0.0 or rho_value <= 0.0:
        r['status'] = 'Both pressure and density must have positive values.'
        return r

    T_initial = None
    results = IAPWS95.temperature_at_fixed_pressure_and_density(P_value, rho_value, T_initial)
    if results is None:
        r['status'] = 'Cannot find a solution, an initial temperature guess may help.'
        return r

    T_value, P_value, rho_value = results

    results = IAPWS95.saturation_pressure_at_fixed_temperature(T_value)
    if results is None:
        kw['Temperature'] = (T_value, 'K')
        kw['Density'] = (rho_value, 'kg/m3')
        return calculate_phase(OVERALL, kw, r)

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rho_value, 'kg/m3')
        
    P_value_sat = results[1]
    rhoV_value = results[2]
    rhoL_value = results[3]
    
    if P_value < P_value_sat:
        r['Quality'][VAPOR] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(VAPOR, kw, r)
    
    elif P_value > P_value_sat:
        r['Quality'][LIQUID] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(LIQUID, kw, r)
    else:
        kw['Density'] = (rhoV_value, 'kg/m3') 
        r = calculate_phase(VAPOR, kw, r)

        kw['Density'] = (rhoL_value, 'kg/m3')
        r = calculate_phase(LIQUID, kw, r)
        r['status'] = 'Results available'        
        return r


def calculate_Pressure_Temperature(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2
    
    r = new_results()
    P_value, P_unit = kw['Pressure']
    T_value, T_unit = kw['Temperature']
    
    if str(P_unit) != 'kPa' or str(T_unit) != 'K':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if T_value < 0.0 or P_value < 0.0:
        r['status'] = 'Both temperature and pressure must be positive.'
        return r

    rho_initial = None
    results = IAPWS95.density_at_fixed_pressure_and_temperature(T_value, P_value, rho_initial)
    if results is None:
        r['status'] = 'Cannot find a solution, an initial density guess may help.'
        return r

    rho_value = results[2]

    results = IAPWS95.saturation_pressure_at_fixed_temperature(T_value)
    if results is None:
        kw['Temperature'] = (T_value, 'K')
        kw['Density'] = (rho_value, 'kg/m3')
        return calculate_phase(OVERALL, kw, r)

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rho_value, 'kg/m3')
        
    P_value_sat = results[1]
    rhoV_value = results[2]
    rhoL_value = results[3]
    
    if P_value < P_value_sat:
        r['Quality'][VAPOR] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(VAPOR, kw, r)
    
    elif P_value > P_value_sat:
        r['Quality'][LIQUID] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(LIQUID, kw, r)
    else:
        kw['Density'] = (rhoV_value, 'kg/m3') 
        r = calculate_phase(VAPOR, kw, r)

        kw['Density'] = (rhoL_value, 'kg/m3')
        r = calculate_phase(LIQUID, kw, r)
        r['status'] = 'Results available'        
        return r
    
    r['status'] = 'Results available'
    return r

def calculate_Enthalpy_Pressure(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2
    
    r = new_results()
    P_value, P_unit = kw['Pressure']
    h_value, h_unit = kw['Enthalpy']
    
    if str(P_unit) != 'kPa' or str(h_unit) != 'kJ/kg':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if P_value < 0.0:
        r['status'] = 'Pressure must be positive.'
        return r

    results = IAPWS95.isenthalpic_at_fixed_pressure(P_value, h_value)
    if results is None:
        r['status'] = 'Cannot find a solution, initial density and temperature guesses may help.'
        return r

    T_value, P_value, rho_value, h_value = results
    
    results = IAPWS95.saturation_pressure_at_fixed_temperature(T_value)
    if results is None:
        kw['Temperature'] = (T_value, 'K')
        kw['Density'] = (rho_value, 'kg/m3')
        return calculate_phase(OVERALL, kw, r)
    
    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rho_value, 'kg/m3')
        
    P_value_sat = results[1]
    rhoV_value = results[2]
    rhoL_value = results[3]
    
    if P_value < P_value_sat:
        r['Quality'][VAPOR] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(VAPOR, kw, r)
    
    elif P_value > P_value_sat:
        r['Quality'][LIQUID] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(LIQUID, kw, r)
    else:
        r['status'] = 'Likely two-phase result, use a quality-enthalpy flash instead'        
        return r
    
    r['status'] = 'Results available'
    return r


def calculate_Enthalpy_Temperature(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2
    
    r = new_results()
    T_value, T_unit = kw['Temperature']
    h_value, h_unit = kw['Enthalpy']
    
    if str(T_unit) != 'K' or str(h_unit) != 'kJ/kg':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if T_value < 0.0:
        r['status'] = 'Temperature must be positive.'
        return r

    results = IAPWS95.isenthalpic_at_fixed_temperature(T_value, h_value)
    if results is None:
        r['status'] = 'Cannot find a solution, initial density and temperature guesses may help.'
        return r

    T_value, P_value, rho_value, h_value = results
    
    results = IAPWS95.saturation_pressure_at_fixed_temperature(T_value)
    if results is None:
        kw['Temperature'] = (T_value, 'K')
        kw['Density'] = (rho_value, 'kg/m3')
        return calculate_phase(OVERALL, kw, r)
    
    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rho_value, 'kg/m3')
        
    P_value_sat = results[1]
    rhoV_value = results[2]
    rhoL_value = results[3]
    
    if P_value < P_value_sat:
        r['Quality'][VAPOR] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(VAPOR, kw, r)
    
    elif P_value > P_value_sat:
        r['Quality'][LIQUID] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(LIQUID, kw, r)
    else:
        r['status'] = 'Likely two-phase result, use a quality-enthalpy flash instead'        
        return r
    
    r['status'] = 'Results available'
    return r


def calculate_Enthalpy_Quality(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2
    
    r = new_results()
    x_value, x_unit = kw['Quality']
    h_value, h_unit = kw['Enthalpy']
    
    if str(x_unit) != 'Unitless' or str(h_unit) != 'kJ/kg':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if x_value < 0.0 or x_value > 1.0:
        r['status'] = 'Quality must be between 0 and 1'
        return r

    T_initial = None
    try:
        results = IAPWS95.isenthalpic_at_fixed_quality(h_value, x_value, T_initial)
    except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError):
        results = None
        
    if results is None:
        r['status'] = 'Cannot find a solution, physical solution unlikely.'
        return r

    T_value, P_value, rhoV_value, rhoL_value, hV_value, hL_value = results
    
    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoV_value, 'kg/m3')
    r['Quality'][VAPOR] = x_value
    r = calculate_phase(VAPOR, kw, r)

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoL_value, 'kg/m3')
    r['Quality'][LIQUID] = 1.0 - x_value
    r = calculate_phase(LIQUID, kw, r)
    r = calculate_overall_phase(r)
        
    r['status'] = 'Results available'
    return r



def calculate_Entropy_Pressure(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2
    
    r = new_results()
    P_value, P_unit = kw['Pressure']
    s_value, s_unit = kw['Entropy']
    
    if str(P_unit) != 'kPa' or str(s_unit) != 'kJ/kg-K':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if P_value < 0.0:
        r['status'] = 'Pressure must be positive.'
        return r

    T_initial = None
    results = IAPWS95.isentropic_at_fixed_pressure(P_value, s_value, T_initial)
    if results is None:
        r['status'] = 'Cannot find a solution, initial density and temperature guesses may help.'
        return r

    T_value, P_value, rho_value, s_value = results
    
    results = IAPWS95.saturation_pressure_at_fixed_temperature(T_value)
    if results is None:
        kw['Temperature'] = (T_value, 'K')
        kw['Density'] = (rho_value, 'kg/m3')
        r['status'] = 'Results available'
        return calculate_phase(OVERALL, kw, r)
    
    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rho_value, 'kg/m3')
        
    P_value_sat = results[1]
    rhoV_value = results[2]
    rhoL_value = results[3]
    
    if P_value < P_value_sat:
        r['Quality'][VAPOR] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(VAPOR, kw, r)
    
    elif P_value > P_value_sat:
        r['Quality'][LIQUID] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(LIQUID, kw, r)
    else:
        r['status'] = 'Likely two-phase result, use a quality-entropy flash instead'        
        return r
    
    r['status'] = 'Results available'
    return r


def calculate_Entropy_Temperature(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2
    
    r = new_results()
    T_value, T_unit = kw['Temperature']
    s_value, s_unit = kw['Entropy']
    
    if str(T_unit) != 'K' or str(s_unit) != 'kJ/kg-K':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if T_value < 0.0:
        r['status'] = 'Temperature must be positive.'
        return r

    results = IAPWS95.isentropic_at_fixed_temperature(T_value, s_value)
    if results is None:
        r['status'] = 'Cannot find a solution, initial density and temperature guesses may help.'
        return r

    T_value, P_value, rho_value, s_value = results
    
    results = IAPWS95.saturation_pressure_at_fixed_temperature(T_value)
    if results is None:
        kw['Temperature'] = (T_value, 'K')
        kw['Density'] = (rho_value, 'kg/m3')
        return calculate_phase(OVERALL, kw, r)
    
    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rho_value, 'kg/m3')
        
    P_value_sat = results[1]
    rhoV_value = results[2]
    rhoL_value = results[3]
    
    if P_value < P_value_sat:
        r['Quality'][VAPOR] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(VAPOR, kw, r)
    
    elif P_value > P_value_sat:
        r['Quality'][LIQUID] = 1.0
        r['status'] = 'Results available'
        return calculate_phase(LIQUID, kw, r)
    else:
        r['status'] = 'Likely two-phase result, use a quality-enthalpy flash instead'        
        return r
    
    r['status'] = 'Results available'
    return r



def calculate_Entropy_Quality(kw):
    OVERALL = 0
    VAPOR = 1
    LIQUID = 2
    
    r = new_results()
    x_value, x_unit = kw['Quality']
    s_value, s_unit = kw['Entropy']
    
    if str(x_unit) != 'Unitless' or str(s_unit) != 'kJ/kg-K':
        r['status'] = 'Internal inconsistency, please report this error.'
        return r

    if x_value < 0.0 or x_value > 1.0:
        r['status'] = 'Quality must be between 0 and 1'
        return r

    T_initial = None
    try:
        results = IAPWS95.isentropic_at_fixed_quality(s_value, x_value, T_initial)
    except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError):
        results = None

    if results is None:
        r['status'] = 'Cannot find a solution, physical solution unlikely.'
        return r

    T_value, P_value, rhoV_value, rhoL_value, sV_value, sL_value = results
    
    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoV_value, 'kg/m3')
    r['Quality'][VAPOR] = x_value
    r = calculate_phase(VAPOR, kw, r)

    kw['Temperature'] = (T_value, 'K')
    kw['Density'] = (rhoL_value, 'kg/m3')
    r['Quality'][LIQUID] = 1.0 - x_value
    r = calculate_phase(LIQUID, kw, r)
    r = calculate_overall_phase(r)
        
    r['status'] = 'Results available'
    return r




def main():
    application = webapp.WSGIApplication([('/', MainHandler),
                                          ('/calc', CalcHandler)], debug=True)
    
    util.run_wsgi_app(application)


if __name__ == '__main__':
    main()

