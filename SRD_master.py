#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:19:00 2023

This function does all the required calculations to map a real fluid to an 
 MPCD fluid. 

@author: lukedebono
"""
#%%
import os
import numpy as np

import matplotlib.pyplot as plt
import regex as re
import pandas as pd

plt.rcParams.update(plt.rcParamsDefault)
#plt.rcParams['text.usetex'] = True
from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
#import seaborn as sns
import math as m
import scipy.stats
from datetime import datetime
from units_lj_scalings import *
from SRD_a_minus_gamma_funcs import *


import Sc_num_est
from box_size_dim_and_integer_bin_constraint import box_size_dim_and_integer_SRD_bin_count_constraint_func
from SRD_calculations_mod import *

# from "module" import "function"

def SRD_MASTER_calc_(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,scaled_timestep,rtol,nu_s,Solvent_bead_SRD_box_density_cp_1, box_size_vec,box_size_vec_nd,SRD_box_size_wrt_solid_beads_check,box_side_length_scaled,T_K,SRD_mass_scale_parameter,lengthscale_parameter,k_b,rho_s,eta_s):
    
    
    
    ### non-dimensionalising inputs 
    scalings_calculation= units_lj_scalings_(SRD_mass_scale_parameter,lengthscale_parameter,k_b,rho_s,eta_s,T_K)

    #energy_parameter=scalings_calculation[0]
    timescale_parameter=scalings_calculation[1]
    temperature_parameter=scalings_calculation[2]
    #scaled_dynamic_viscosity=scalings_calculation[3]
    scaled_nu_s=(scalings_calculation[4])
    #scaled_rho_s=scalings_calculation[5]
    scaled_temp=T_K/temperature_parameter

#########################################################################
    ## box size calculations 
    
    #box_size_vec = np.array([box_side_length/number_boxes_vec])
    #box_size_vec_nd=np.array([box_side_length_scaled/number_boxes_vec])
    SRD_box_size_wrt_solid_beads =box_size_vec_nd
    SRD_box_size_wrt_solid_beads_check = box_size_vec
#########################################################################    
    
    ## Calculating gamma values
    gammas=SRD_a_minus_gamma_funcs_(Solvent_bead_SRD_box_density_cp_1)
    gamma_1=gammas[0]
    gamma_2=gammas[1]
#########################################################################   


    
    SRD_timesteps =SRD_timestep_non_dimensional_and_dimensional(box_side_length,number_boxes_vec,box_side_length_scaled,Solvent_bead_SRD_box_density_cp_1,mass_fluid_particle_wrt_pf_cp_mthd_1,SRD_mass_scale_parameter,nu_s,k_b,scaled_nu_s,scaled_temp,lengthscale_parameter,timescale_parameter,T_K,gamma_1,gamma_2)
    SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_timesteps[0]
    SRD_timestep_cp_1_based_on_sphere_pf_neg_nd=SRD_timesteps[1]
    SRD_timestep_cp_1_based_on_sphere_pf_neg_re_dim=SRD_timesteps[2]
    SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_timesteps[3]
    SRD_timestep_cp_1_based_on_sphere_pf_pos_nd=SRD_timesteps[4]
    SRD_timestep_cp_1_based_on_sphere_pf_pos_re_dim=SRD_timesteps[5]
#########################################################################
    #this section checks if the nd calc matches the dim calc 
    arr=SRD_timestep_cp_1_based_on_sphere_pf_neg_re_dim
    arr1=SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check
    print("Negative n-d timestep solution matches d solution is T/F?")
    print((np.allclose(arr,arr1,atol=0.000001, equal_nan=True)))
    if np.allclose(arr,arr1,atol=0.000001, equal_nan=True) == 1:
        print("Neg solution calculation Success") 
    else: 
        print("Neg solution calculation failure")
        breakpoint()
    
       
    arr=SRD_timestep_cp_1_based_on_sphere_pf_pos_re_dim
    arr1=SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check
    print("Positve n-d timestep solution matches d solution is T/F?")
    print((np.allclose(arr,arr1,atol=0.000001, equal_nan=True)))
    if np.allclose(arr,arr1,atol=0.000001, equal_nan=True) == 1:
        print("Pos solution calculation Success") 
    else: 
        print("Pos solution calculation failure")
        breakpoint()
      
#########################################################################    
    
    ## schmidt number calculation from dimensional numbers
    sc_pos_soln= Sc_num_est.Schmidt_num_est(gamma_1, gamma_2, k_b, T_K, SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check, SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check, rho_s, Solvent_bead_SRD_box_density_cp_1, SRD_box_size_wrt_solid_beads_check)[0]       
    sc_neg_soln= Sc_num_est.Schmidt_num_est(gamma_1, gamma_2, k_b, T_K, SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check, SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check, rho_s, Solvent_bead_SRD_box_density_cp_1, SRD_box_size_wrt_solid_beads_check)[1]       
#########################################################################           
    ## calculating mean free path 
    
    mean_free_paths= SRD_mean_free_path(mass_fluid_particle_wrt_pf_cp_mthd_1,SRD_mass_scale_parameter,scaled_temp,k_b,T_K,SRD_timestep_cp_1_based_on_sphere_pf_pos_nd,SRD_timestep_cp_1_based_on_sphere_pf_neg_nd,SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check,SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check,lengthscale_parameter)
    mean_free_path_pf_SRD_particles_cp_mthd_1_neg_dim_check=mean_free_paths[0]
    mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=mean_free_paths[1]
    mean_free_path_pf_SRD_particles_cp_mthd_1_neg_re_dim=mean_free_paths[2]
    mean_free_path_pf_SRD_particles_cp_mthd_1_pos_dim_check=mean_free_paths[3]
    mean_free_path_pf_SRD_particles_cp_mthd_1_pos_nd=mean_free_paths[4]
    mean_free_path_pf_SRD_particles_cp_mthd_1_pos_re_dim=mean_free_paths[5]
    
    arr=mean_free_path_pf_SRD_particles_cp_mthd_1_neg_re_dim
    arr1=mean_free_path_pf_SRD_particles_cp_mthd_1_neg_dim_check
    print("Negative n-d mfp solution matches d solution is T/F?")
    print((np.allclose(arr,arr1,atol=0.000001, equal_nan=True)))
    if np.allclose(arr,arr1,atol=0.000001, equal_nan=True) == 1:
        print("Neg solution calculation Success") 
    else: 
        print("Neg solution calculation failure")
        breakpoint()
    
       
    arr=mean_free_path_pf_SRD_particles_cp_mthd_1_pos_re_dim

    arr1=mean_free_path_pf_SRD_particles_cp_mthd_1_pos_dim_check
    print("Positve n-d mfp solution matches d solution is T/F?")
    print((np.allclose(arr,arr1,rtol=rtol, equal_nan=True)))
    if np.allclose(arr,arr1,rtol=rtol, equal_nan=True) == 1:
        print("Pos solution calculation Success") 
    else: 
        print("Pos solution calculation failure")
        breakpoint()
    
   
   
   
   
#########################################################################       
    # print('box_side_length_scaled:',box_side_length_scaled) 
    # print('SRD_box_size_wrt_solid_beads:',SRD_box_size_wrt_solid_beads[0,54])   
#     Number_of_SRD_boxes_in_sim_box_wrt_pf = (box_side_length_scaled**3)/(SRD_box_size_wrt_solid_beads**3)  
#    # print('Number_of_SRD_boxes_in_sim_box_wrt_pf:',Number_of_SRD_boxes_in_sim_box_wrt_pf[0,54])
#     number_SRD_particles_wrt_pf_cp_mthd_1_pos = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf[z,:] * Solvent_bead_SRD_box_density_cp_1)
#     number_SRD_particles_wrt_pf_cp_mthd_1_neg = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf[z,:] * Solvent_bead_SRD_box_density_cp_1)
    #print('number_SRD_particles_wrt_pf_cp_mthd_1_neg:',number_SRD_particles_wrt_pf_cp_mthd_1_neg[0,54] )

    Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos = SRD_timestep_cp_1_based_on_sphere_pf_pos_nd/scaled_timestep
    Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg = SRD_timestep_cp_1_based_on_sphere_pf_neg_nd/scaled_timestep

    ## pre-calc to remove timestep ratios which arent close to integer values,
    # reduces error in intial temperature
    check=np.zeros((Solvent_bead_SRD_box_density_cp_1.size,box_size_vec.size))
    check= check.astype('float64')
    comparison_pos = (np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos)-Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos)/Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos
    comparison_neg = (np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)-Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)/Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg
#########################################################################
    
    return  sc_pos_soln,sc_neg_soln,mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd,mean_free_path_pf_SRD_particles_cp_mthd_1_pos_nd,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos,comparison_neg,comparison_pos,SRD_timestep_cp_1_based_on_sphere_pf_neg_nd,SRD_timestep_cp_1_based_on_sphere_pf_pos_nd
 
 
    
