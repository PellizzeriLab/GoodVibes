# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

import ctypes, math, os.path, sys
import numpy as np

# Importing regardless of relative import
try:
    from .io import *
except:
    from io import *

# PHYSICAL CONSTANTS                                      UNITS
GAS_CONSTANT = 8.3144621  # J / K / mol
PLANCK_CONSTANT = 6.62606957e-34  # J * s
BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K
SPEED_OF_LIGHT = 2.99792458e10  # cm / s
AVOGADRO_CONSTANT = 6.0221415e23  # 1 / mol
AMU_to_KG = 1.66053886E-27  # UNIT CONVERSION
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION


def calc_heat_capacity_vib(frequency_wn, temperature, freq_scale_factor):
    """
    Rigid rotor harmonic oscillator (RRHO) heat capacity evaluation - this is the default treatment

    Entropic contributions (J/(mol*K)) according to a rigid-rotor
    harmonic-oscillator description for a list of vibrational modes
    C_v_vib = Sum(R*(hv/kT)^2*((e^(hv/kT)/(e^(hv/kT)-1)^2))
                           
    Parameters:
    frequency_wn (list): list of frequencies parsed from file.
    temperature (float): temperature for calculations to be performed at.
    freq_scale_factor (float): frequency scaling factor based on level of theory and basis set used.
    
    Returns:
    float: RRHO C_V of chemical system.
    """
    
    factor = [(PLANCK_CONSTANT * freq * SPEED_OF_LIGHT * freq_scale_factor) / (BOLTZMANN_CONSTANT * temperature)
                  for freq in frequency_wn]

    C_V_Freq = np.array([GAS_CONSTANT * entry**2 *(math.exp(entry)/(math.exp(entry)-1)**2)
               for entry in factor])


    return C_V_Freq



def calc_heat_capacity_trans():
     """
     Rigid rotor harmonic oscillator (RRHO) heat capacity evaluation - this is the default treatment

     C_v_trans = (3/2)* R
                            
     Parameters:
  
     
     Returns:
     float: translationsal C_V of chemical system.
     """
     
     C_V_Trans = (3/2)*float(GAS_CONSTANT)

     return C_V_Trans
 
def calc_heat_capacity_rot(zpe, linear):
    """
    Rotational energy evaluation

    Calculates the rotational energy (J/mol)
    Cv_rot = 0 (atomic) ; RT (linear); 3/2 RT (non-linear)

    Parameters:

    zpe (float): zero point energy of chemical system.
    linear (bool): flag for linear molecules, changes how calculation is performed.

    Returns:
    float: roational C_V of chemical system.
    """
    if zpe == 0.0:
        C_V_ROT = 0.0
    elif linear == 1:
        C_V_ROT = GAS_CONSTANT
    else:
        C_V_ROT = 1.5 * GAS_CONSTANT
    
    return C_V_ROT

    