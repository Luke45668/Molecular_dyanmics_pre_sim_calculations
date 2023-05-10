#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 13:20:52 2023

This script will do the SRD calculations in dimensional form, then do exactly the same 
calculations in dimensionless form. This will allow me to figure out where the errors are 
coming from 



@author: lukedebono
"""
#%%

import os
import numpy as np

import matplotlib.pyplot as plt
import regex as re
import pandas as pd

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
#import seaborn as sns
import math as m
import scipy.stats
from datetime import datetime

#%% fixed values 
fluid_name='Nitrogen'
tolerance=0.01 # for solution error 
atol=0.01
rtol=0.00001
number_of_test_points =50
Solvent_bead_SRD_box_density_cp_1 = np.array([(np.linspace(3.5,1000,number_of_test_points))])
number_of_M_cp_1=Solvent_bead_SRD_box_density_cp_1.shape[1]
scaled_timestep=0.01

#%% fixed values for nitrogen 
rho_s = 847 #kg/m^3
T_K=72.2 #+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1

# linear interpolation
# to determine true visc from NIST data
eta_s_1 = 0.00022081
rho_1=	843.96
eta_s_2 = 0.00023138
rho_2 = 850.73	
delta_rho=rho_2-rho_1 
delta_eta_s =eta_s_2 - eta_s_1 
grad_eta_s_rho_s = delta_eta_s/delta_rho

eta_s_NIST=0.00022081 + ((rho_s -rho_1)*grad_eta_s_rho_s) 
eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s

#determine side length of simulaton box 
r_particle =50e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)

number_boxes_vec=np.linspace(2,32,31)
box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1.T

#%%
length_multiplier=0.1
lengthscale_parameter = length_multiplier*r_particle
box_side_length_scaled=box_side_length/lengthscale_parameter
box_size_to_lengthscale=box_size_vec/lengthscale_parameter
mass_multiplier=1
SRD_mass_scale_parameter = mass_multiplier* rho_s * (lengthscale_parameter**3)
r_particle_scaled = r_particle/lengthscale_parameter

import units_lj_scalings
scalings_calculation= units_lj_scalings.units_lj_scalings_(SRD_mass_scale_parameter,lengthscale_parameter,k_b,rho_s,eta_s,T_K)

energy_parameter=scalings_calculation[0]
timescale_parameter=scalings_calculation[1]
temperature_parameter=scalings_calculation[2]
scaled_dynamic_viscosity=scalings_calculation[3]
scaled_nu_s=(scalings_calculation[4])
scaled_rho_s=scalings_calculation[5]
scaled_temp=T_K/temperature_parameter

#######################################################################################################
    ### START OF  DIMENSIONAL COMPARISON ###
    ### COMPARE EACH STAGE STEP BY STEP TO SEE WHERE THE ISSUE HAPPENS ###
#######################################################################################################

#%%
# import SRD_master_calc_dimensional
# from units_lj_scalings import *
# from SRD_a_minus_gamma_funcs import *
# import srd_timestep_non_dim
# from  srd_timestep_dim import *
# import Sc_num_est
# from box_size_dim_and_integer_bin_constraint import box_size_dim_and_integer_SRD_bin_count_constraint_func
# from srd_timestep_dim import srd_timestep_dimensional_


# SRD_dimensional_master_data=SRD_master_calc_dimensional.SRD_MASTER_calc_dimensional(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec,T_K,k_b,rho_s,eta_s)
# SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_dimensional_master_data[0]
# SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_dimensional_master_data[1]

# SRD_step_pos_nd_check=SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check/timescale_parameter
# SRD_step_neg_nd_check=SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check/timescale_parameter

# SRD_MD_ratio_neg=SRD_step_neg_nd_check/scaled_timestep
# SRD_MD_ratio_pos=SRD_step_pos_nd_check/scaled_timestep
# sc_pos_soln_from_dim=SRD_dimensional_master_data[2] #yes 
# sc_neg_soln_from_dim=SRD_dimensional_master_data[3] #yes 


# mean_free_path_pf_SRD_particles_cp_mthd_1_neg_from_dim=(SRD_dimensional_master_data[4])
# mean_free_path_pf_SRD_particles_cp_mthd_1_pos_from_dim=(SRD_dimensional_master_data[5])
# mean_free_path_to_box_ratio_neg_from_dim=mean_free_path_pf_SRD_particles_cp_mthd_1_neg_from_dim/box_size_vec
# mean_free_path_to_box_ratio_pos_from_dim=mean_free_path_pf_SRD_particles_cp_mthd_1_pos_from_dim/box_size_vec

# number_SRD_particles_wrt_pf_cp_mthd_1_neg_from_dim=SRD_dimensional_master_data[6]
# number_SRD_particles_wrt_pf_cp_mthd_1_pos_from_dim=SRD_dimensional_master_data[7]
# mass_fluid_particle_wrt_pf_cp_mthd_1_from_dim=SRD_dimensional_master_data[8]

#######################################################################################################
    ### START OF NON DIMENSIONAL COMPARISON ###
    ### COMPARE EACH STAGE STEP BY STEP TO SEE WHERE THE ISSUE HAPPENS ###
#######################################################################################################

#%%  
import numpy as np
from SRD_master import *

box_size_vec = np.array([box_side_length/number_boxes_vec])
box_size_vec_nd=np.array([box_side_length_scaled/number_boxes_vec])
number_of_boxes_in_each_dim=number_boxes_vec

#gamma_1 = ((1 - ((1-np.exp(-(Solvent_bead_SRD_box_density_cp_1)))/(Solvent_bead_SRD_box_density_cp_1))).T)
#gamma_2 = ( (((Solvent_bead_SRD_box_density_cp_1)+2)/((Solvent_bead_SRD_box_density_cp_1)-1)).T)


SRD_non_dimensional_master_data=SRD_MASTER_calc_(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,scaled_timestep,atol,rtol,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,box_side_length_scaled,T_K,SRD_mass_scale_parameter,lengthscale_parameter,energy_parameter,k_b,rho_s,eta_s)
sc_pos_soln=SRD_non_dimensional_master_data[0]
sc_neg_soln=SRD_non_dimensional_master_data[1]

mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=SRD_non_dimensional_master_data[2]
mean_free_path_pf_SRD_particles_cp_mthd_1_pos_nd=SRD_non_dimensional_master_data[3]

Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg=SRD_non_dimensional_master_data[4]
Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos=SRD_non_dimensional_master_data[5]

number_SRD_particles_wrt_pf_cp_mthd_1_neg=SRD_non_dimensional_master_data[6]
number_SRD_particles_wrt_pf_cp_mthd_1_pos=SRD_non_dimensional_master_data[7]

mass_fluid_particle_wrt_pf_cp_mthd_1=SRD_non_dimensional_master_data[8]

comparison_neg=SRD_non_dimensional_master_data[9]
comparison_pos=SRD_non_dimensional_master_data[10]

SRD_timestep_cp_1_based_on_sphere_pf_neg_nd=SRD_non_dimensional_master_data[11]
SRD_timestep_cp_1_based_on_sphere_pf_pos_nd=SRD_non_dimensional_master_data[12]
 
 
    


#%%
number_boxes_vec=np.linspace(2,32,31)
box_size_vec = np.array([box_side_length/number_boxes_vec])
box_size_vec_nd=np.array([box_side_length_scaled/number_boxes_vec])
number_of_boxes_in_each_dim=number_boxes_vec
SRD_box_size_wrt_solid_beads =box_size_vec_nd #np.array([box_side_length_scaled/number_of_boxes_in_each_dim])
gamma_1_1= ((1 - ((1-np.exp(-(Solvent_bead_SRD_box_density_cp_1)))/(Solvent_bead_SRD_box_density_cp_1))).T)
gamma_2_2 = ( (((Solvent_bead_SRD_box_density_cp_1)+2)/((Solvent_bead_SRD_box_density_cp_1)-1)).T)
gammas=SRD_a_minus_gamma_funcs_(Solvent_bead_SRD_box_density_cp_1)
gamma_1=gammas[0]
gamma_2=gammas[1]
SRD_box_size_wrt_solid_beads_check = box_size_vec

# %%
# splitting the engative calculation in the exact same way
#dim 
srd1_neg =nu_s - np.sqrt((nu_s**2) - ((k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1)))
srd1a_neg= np.sqrt((nu_s**2) - ((k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1)))


srd2=(k_b*T_K*gamma_2) /(2*mass_fluid_particle_wrt_pf_cp_mthd_1)

out_SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check = (srd1_neg/srd2)

#non-dim 
mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/SRD_mass_scale_parameter


srd1nd_neg=scaled_nu_s-np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))
srd1and_neg =np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(lengthscale_parameter*SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))


srd2nd = (scaled_temp * gamma_2 )/(2*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)
out_SRD_timestep_cp_1_based_on_sphere_pf_neg_nd=srd1nd_neg/srd2nd
out_SRD_timestep_cp_1_based_on_sphere_pf_neg_re_dim=out_SRD_timestep_cp_1_based_on_sphere_pf_neg_nd*timescale_parameter
# %%# %%
# splitting the positive calculation in the exact same way
#dim 
srd1_pos =nu_s + np.sqrt((nu_s**2) - ((k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1)))
srd1a_pos=np.sqrt((nu_s**2) - ((k_b*T_K*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads_check**2))/(18*mass_fluid_particle_wrt_pf_cp_mthd_1)))




out_SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check = (srd1_pos/srd2)

#non-dim 
mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/SRD_mass_scale_parameter


srd1nd_pos=scaled_nu_s+np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))
srd1and_pos=np.sqrt((scaled_nu_s**2)-((scaled_temp*gamma_1*gamma_2*(lengthscale_parameter*SRD_box_size_wrt_solid_beads)**2)/(18*mass_fluid_particle_wrt_pf_cp_mthd_1_nd)))



out_SRD_timestep_cp_1_based_on_sphere_pf_pos_nd=srd1nd_pos/srd2nd
out_SRD_timestep_cp_1_based_on_sphere_pf_pos_re_dim=out_SRD_timestep_cp_1_based_on_sphere_pf_pos_nd*timescale_parameter
# %%
