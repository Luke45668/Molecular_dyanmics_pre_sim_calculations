#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:34:33 2023

@author: lukedebono
This function calculates the srd timestep from non-dimensional parameters

"""
import numpy as np 
def srd_timestep_non_dim_(mass_fluid_particle_wrt_pf_cp_mthd_1,scaled_nu_s,SRD_mass_scale_parameter,lengthscale_parameter,scaled_temp,gamma_1,gamma_2,SRD_box_size_wrt_solid_beads):  
    mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/SRD_mass_scale_parameter
    srd1nd_neg=scaled_nu_s-np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))
    srd1nd_pos=scaled_nu_s+np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))
    srd2nd = (scaled_temp * gamma_2 )/(2*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)
    SRD_timestep_cp_1_based_on_sphere_pf_neg_nd=srd1nd_neg/srd2nd
    SRD_timestep_cp_1_based_on_sphere_pf_pos_nd=srd1nd_pos/srd2nd
    
    
    return SRD_timestep_cp_1_based_on_sphere_pf_neg_nd,SRD_timestep_cp_1_based_on_sphere_pf_pos_nd