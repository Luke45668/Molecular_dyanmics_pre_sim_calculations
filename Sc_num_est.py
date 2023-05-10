#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:22:19 2023

@author: lukedebono
"""

import numpy as np 


def Schmidt_num_est(gamma_1,gamma_2,k_b,T_K,SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check,SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check,rho_s,Solvent_bead_SRD_box_density_cp_1,SRD_box_size_wrt_solid_beads_check):
    M=Solvent_bead_SRD_box_density_cp_1 
    gamma_3 = ((3*M /(M-1 + np.exp(-M))) -1 )
    D_f_neg_dm =( M*k_b*T_K *SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check.T * gamma_3)/ (rho_s * SRD_box_size_wrt_solid_beads_check**3).T 
    D_f_pos_dm=( M*k_b*T_K *SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check.T * gamma_3)/ (rho_s * SRD_box_size_wrt_solid_beads_check**3).T
    nu_neg_dm = (SRD_box_size_wrt_solid_beads_check**2/18*SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check) * gamma_1 + (k_b*T_K*SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check*gamma_2)/(4*rho_s * (SRD_box_size_wrt_solid_beads_check**3))
    nu_pos_dm = (SRD_box_size_wrt_solid_beads_check**2/18*SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check) * gamma_1 + (k_b*T_K*SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check*gamma_2)/(4*rho_s * (SRD_box_size_wrt_solid_beads_check**3))
    
    Sc_pos=nu_pos_dm/D_f_pos_dm.T
    Sc_neg=nu_neg_dm/D_f_neg_dm.T
    
    return Sc_pos,Sc_neg,D_f_neg_dm,D_f_pos_dm

    
    
   