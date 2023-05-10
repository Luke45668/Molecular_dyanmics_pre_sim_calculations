#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 23 2023
This script does both non dimension and dimensional SRD timesteo calculations 
@author: lukedebono
"""
import os
import numpy as np





def SRD_timestep_non_dimensional_and_dimensional(box_side_length,number_boxes_vec,box_side_length_scaled,Solvent_bead_SRD_box_density_cp_1,mass_fluid_particle_wrt_pf_cp_mthd_1,SRD_mass_scale_parameter,nu_s,k_b,scaled_nu_s,scaled_temp,lengthscale_parameter,timescale_parameter,T_K,gamma_1,gamma_2):
    #print(box_side_length,number_boxes_vec,box_side_length_scaled,Solvent_bead_SRD_box_density_cp_1,mass_fluid_particle_wrt_pf_cp_mthd_1,SRD_mass_scale_parameter,nu_s,k_b,scaled_nu_s,scaled_temp,lengthscale_parameter,timescale_parameter,T_K,gamma_1,gamma_2)
    
    box_size_vec = np.array([box_side_length/number_boxes_vec])
    box_size_vec_nd=np.array([box_side_length_scaled/number_boxes_vec])
    number_of_boxes_in_each_dim=number_boxes_vec
    SRD_box_size_wrt_solid_beads =box_size_vec_nd 
    #gamma_1 = ((1 - ((1-np.exp(-(Solvent_bead_SRD_box_density_cp_1)))/(Solvent_bead_SRD_box_density_cp_1))).T)
    #gamma_2 = ( (((Solvent_bead_SRD_box_density_cp_1)+2)/((Solvent_bead_SRD_box_density_cp_1)-1)).T)
    SRD_box_size_wrt_solid_beads_check = box_size_vec
    srd2=(k_b*T_K*gamma_2) /(2*mass_fluid_particle_wrt_pf_cp_mthd_1)
    
    #### negative dimensional 
   
    srd1_neg =nu_s - np.sqrt((nu_s**2) - ((k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1)))

    SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check = (srd1_neg/srd2)
    
    #non-dim 
    mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/SRD_mass_scale_parameter
    
    srd1nd_neg=scaled_nu_s-np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))
    
    srd2nd = (scaled_temp * gamma_2 )/(2*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)
    
    SRD_timestep_cp_1_based_on_sphere_pf_neg_nd=srd1nd_neg/srd2nd
    
    # negative re-dim check 
    SRD_timestep_cp_1_based_on_sphere_pf_neg_re_dim=SRD_timestep_cp_1_based_on_sphere_pf_neg_nd*timescale_parameter
    
    
    #dim pos  
    srd1_pos =nu_s + np.sqrt((nu_s**2) - ((k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1)))
    SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check = (srd1_pos/srd2)

    #non-dim pos 

    srd1nd_pos=scaled_nu_s+np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))

    SRD_timestep_cp_1_based_on_sphere_pf_pos_nd=srd1nd_pos/srd2nd
    # pos re-dim check 
    
    SRD_timestep_cp_1_based_on_sphere_pf_pos_re_dim=SRD_timestep_cp_1_based_on_sphere_pf_pos_nd*timescale_parameter
    
    return SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check,SRD_timestep_cp_1_based_on_sphere_pf_neg_nd,SRD_timestep_cp_1_based_on_sphere_pf_neg_re_dim,SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check,SRD_timestep_cp_1_based_on_sphere_pf_pos_nd,SRD_timestep_cp_1_based_on_sphere_pf_pos_re_dim


def SRD_mean_free_path(mass_fluid_particle_wrt_pf_cp_mthd_1,SRD_mass_scale_parameter,scaled_temp,k_b,T_K,SRD_timestep_cp_1_based_on_sphere_pf_pos_nd,SRD_timestep_cp_1_based_on_sphere_pf_neg_nd,SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check,SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check,lengthscale_parameter): 
    mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/SRD_mass_scale_parameter

    mean_free_path_pf_SRD_particles_cp_mthd_1_pos_nd =SRD_timestep_cp_1_based_on_sphere_pf_pos_nd*np.sqrt((scaled_temp)/(mass_fluid_particle_wrt_pf_cp_mthd_1_nd))
    mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd =SRD_timestep_cp_1_based_on_sphere_pf_neg_nd*np.sqrt((scaled_temp)/(mass_fluid_particle_wrt_pf_cp_mthd_1_nd))
    
    mean_free_path_pf_SRD_particles_cp_mthd_1_pos_dim_check =SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check*np.sqrt((T_K*k_b)/(mass_fluid_particle_wrt_pf_cp_mthd_1))
    mean_free_path_pf_SRD_particles_cp_mthd_1_neg_dim_check =SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check*np.sqrt((T_K*k_b)/(mass_fluid_particle_wrt_pf_cp_mthd_1))
    
    mean_free_path_pf_SRD_particles_cp_mthd_1_pos_re_dim= mean_free_path_pf_SRD_particles_cp_mthd_1_pos_nd*lengthscale_parameter
    mean_free_path_pf_SRD_particles_cp_mthd_1_neg_re_dim= mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd*lengthscale_parameter
    
    return mean_free_path_pf_SRD_particles_cp_mthd_1_neg_dim_check,mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd,mean_free_path_pf_SRD_particles_cp_mthd_1_neg_re_dim,mean_free_path_pf_SRD_particles_cp_mthd_1_pos_dim_check,mean_free_path_pf_SRD_particles_cp_mthd_1_pos_nd, mean_free_path_pf_SRD_particles_cp_mthd_1_pos_re_dim
    
    
    
    