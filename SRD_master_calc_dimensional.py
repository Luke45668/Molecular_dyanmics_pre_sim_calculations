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
import SRD_a_minus_gamma_funcs
from srd_timestep_dim import srd_timestep_dimensional_
import Sc_num_est
#%%

def SRD_MASTER_calc_dimensional(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,T_K,k_b,rho_s,eta_s):
    
    
    
#########################################################################
    ## box size calculations 
    # number_boxes_vec=box_side_length_scaled/box_size_vec 
    # number_boxes_vec_round=np.round(number_boxes_vec)
    # box_size_vec=box_side_length_scaled/number_boxes_vec_round
    # box_side_length=box_side_length_scaled *lengthscale_parameter
    # number_boxes_vec=box_side_length_scaled/box_size_vec 
    # r_particle_scaled=r_particle/lengthscale_parameter
    number_of_boxes_in_each_dim=number_boxes_vec
    # SRD_box_size_wrt_solid_beads = box_side_length_scaled/number_of_boxes_in_each_dim
    SRD_box_size_wrt_solid_beads_check = box_size_vec
#########################################################################    
    

    ## Calculating gamma values
    gamma_1=SRD_a_minus_gamma_funcs.SRD_a_minus_gamma_funcs_(Solvent_bead_SRD_box_density_cp_1)[0]
    gamma_2=SRD_a_minus_gamma_funcs.SRD_a_minus_gamma_funcs_(Solvent_bead_SRD_box_density_cp_1)[1]
#########################################################################   
    
    ## Calculatiung mass of SRD particles
    volume_solid_particle = (4/3)*np.pi*(r_particle**3)
    volume_solvent_without_solid = box_side_length**(3)
    #volume_solvent_with_solid = box_side_length_scaled**(3)- volume_solid_particle
    #without solid
    Number_of_SRD_boxes_in_sim_box_wrt_pf = (box_side_length**3)/(SRD_box_size_wrt_solid_beads_check**3)
    number_SRD_particles_wrt_pf_cp_mthd_1 = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf * Solvent_bead_SRD_box_density_cp_1.T)
    #mass_fluid_particle_wrt_pf_cp_mthd_1= (volume_solvent_without_solid*rho_s)/number_SRD_particles_wrt_pf_cp_mthd_1
    #with solid 
    # add this section 

    
    ## calculating SRD timestep dimensional
    SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check =  srd_timestep_dimensional_(mass_fluid_particle_wrt_pf_cp_mthd_1,nu_s,k_b,T_K,gamma_1,gamma_2,rho_s,SRD_box_size_wrt_solid_beads_check)[0]
    SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check =   srd_timestep_dimensional_(mass_fluid_particle_wrt_pf_cp_mthd_1,nu_s,k_b,T_K,gamma_1,gamma_2,rho_s,SRD_box_size_wrt_solid_beads_check)[1]
    
#########################################################################    
    
    ## schmidt number calculation from dimensional numbers
    sc_pos_soln= Sc_num_est.Schmidt_num_est(gamma_1, gamma_2, k_b, T_K, SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check, SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check, rho_s, Solvent_bead_SRD_box_density_cp_1, SRD_box_size_wrt_solid_beads_check)[0]       
    sc_neg_soln= Sc_num_est.Schmidt_num_est(gamma_1, gamma_2, k_b, T_K, SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check, SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check, rho_s, Solvent_bead_SRD_box_density_cp_1, SRD_box_size_wrt_solid_beads_check)[1]       
#########################################################################           
    ## calculating mean free path 
    mean_free_path_pf_SRD_particles_cp_mthd_1_pos =SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check* np.sqrt((k_b*T_K)/(mass_fluid_particle_wrt_pf_cp_mthd_1))
    mean_free_path_pf_SRD_particles_cp_mthd_1_neg = SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check* np.sqrt((k_b*T_K)/(mass_fluid_particle_wrt_pf_cp_mthd_1))

    number_SRD_particles_wrt_pf_cp_mthd_1_pos = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf * Solvent_bead_SRD_box_density_cp_1.T)
    number_SRD_particles_wrt_pf_cp_mthd_1_neg = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf * Solvent_bead_SRD_box_density_cp_1.T)

  
    return  SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check, SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check,sc_pos_soln,sc_neg_soln,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,mean_free_path_pf_SRD_particles_cp_mthd_1_pos,number_SRD_particles_wrt_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_pos,mass_fluid_particle_wrt_pf_cp_mthd_1
 
 
    
# %%
