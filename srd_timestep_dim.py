#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:39:47 2023

@author: lukedebono
This function calculates the dimensional SRD timesteps 
"""
import numpy as np 
def srd_timestep_dimensional_(mass_fluid_particle_wrt_pf_cp_mthd_1,nu_s,k_b,T_K,gamma_1,gamma_2,rho_s,SRD_box_size_wrt_solid_beads_check):
    srd1_neg =nu_s - np.sqrt((nu_s**2) - ((k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1)))
    srd1_pos =nu_s + np.sqrt((nu_s**2) - (k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1))
    srd2=(k_b*T_K*gamma_2) /(2*mass_fluid_particle_wrt_pf_cp_mthd_1)
    
    
    SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check = (srd1_neg/srd2 )
    SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check = (srd1_pos/srd2 )
    return SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check,SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check
