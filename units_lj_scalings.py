#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 09:04:52 2023

@author: lukedebono
This function will take dimensionless scaling inputs and ensure conformation with
the lammps units lj scalings 
"""

import numpy as np

def units_lj_scalings_(SRD_mass_scale_parameter,lengthscale_parameter,k_b,rho_s,eta_s,T_K):
    #energy_parameter = (SRD_mass_scale_parameter * (lengthscale_parameter**2))/(timescale_parameter**2)
    energy_parameter=T_K*k_b
    timescale_parameter= np.sqrt((SRD_mass_scale_parameter*(lengthscale_parameter**2)/energy_parameter))
    temperature_parameter =energy_parameter/k_b
    scaled_rho_s = rho_s*(lengthscale_parameter**3)/SRD_mass_scale_parameter
    scaled_dynamic_viscosity = eta_s * (lengthscale_parameter**(3))/ (energy_parameter* timescale_parameter)
    scaled_nu_s = np.float128( (eta_s/rho_s)*(SRD_mass_scale_parameter/(energy_parameter*timescale_parameter)  ) )#scaled_dynamic_viscosity/scaled_rho_s
    return energy_parameter,timescale_parameter,temperature_parameter,scaled_dynamic_viscosity,scaled_nu_s,scaled_rho_s
