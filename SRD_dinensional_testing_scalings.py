#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 30 2023

This script will do the SRD calculations in dimensional form, which will allow us to pick scalings which give 
order 1 mean free paths, like in Ripoll and Petersen 



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


#%%
fluid_name='Nitrogen'
tolerance=0.01
number_of_test_points =10
Solvent_bead_SRD_box_density_cp_1 = np.array([(np.linspace(10,100,number_of_test_points))])
number_of_M_cp_1=Solvent_bead_SRD_box_density_cp_1.shape[1]
number_boxes_var=100
min_number_boxes_for_particle_size=6
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)

#determine side length of simulaton box 
r_particle =50e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)

### physical parameters for Nitrogen
# fixed values for nitrogen 
rho_s = 847 #kg/m^3
Temp_visc_multiplier= 1#0.0000099 #0.00005
T_K=72.2 *  Temp_visc_multiplier#+273.15 #Kelvin
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
eta_s_multiplier=1
eta_s=eta_s_NIST*Temp_visc_multiplier #*1000 to convert kg to g

nu_s = (eta_s/rho_s) 
box_size_vec = np.array([box_side_length/number_boxes_vec])

mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1.T
#%% dimensional SRD calculations 

from SRD_master_calc_dimensional import * 

SRD_dimensional_data= SRD_MASTER_calc_dimensional(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,T_K,k_b,rho_s,eta_s)
 
SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_dimensional_data[0]
SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_dimensional_data[1]
sc_pos_soln=SRD_dimensional_data[2]
sc_neg_soln=SRD_dimensional_data[3]
mean_free_path_pf_SRD_particles_cp_mthd_1_neg=SRD_dimensional_data[4]
mean_free_path_pf_SRD_particles_cp_mthd_1_pos=SRD_dimensional_data[5]
number_SRD_particles_wrt_pf_cp_mthd_1_neg=SRD_dimensional_data[6]
number_SRD_particles_wrt_pf_cp_mthd_1_pos=SRD_dimensional_data[7]
mass_fluid_particle_wrt_pf_cp_mthd_1=SRD_dimensional_data[8]
 

 
 
# %%
scaled_timestep=0.01
lengthscale=r_particle * 0.01
r_particle_scaled=r_particle/lengthscale

box_size_vec_nd= box_size_vec/lengthscale
mass_scale= 10000000*rho_s* (lengthscale**3)
mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/mass_scale
energy_scale=k_b*T_K
timescale= np.sqrt((mass_scale*(lengthscale**2))/energy_scale)
scaled_temp=T_K/(energy_scale/k_b)
SRD_timestep_cp_1_based_on_sphere_neg_nd= SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check/timescale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=mean_free_path_pf_SRD_particles_cp_mthd_1_neg/lengthscale
mean_free_path_to_box_ratio_neg= mean_free_path_pf_SRD_particles_cp_mthd_1_neg/box_size_vec
mean_free_path_to_box_ratio_pos= mean_free_path_pf_SRD_particles_cp_mthd_1_pos/box_size_vec
SRD_MD_Ratio_neg= SRD_timestep_cp_1_based_on_sphere_neg_nd/scaled_timestep
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
from petersen_plotting import *
plt.rcParams.update({'font.size': 25})
def sc_vs_mfp_to_collision_cell(mean_free_path_to_box_ratio_pos,mean_free_path_to_box_ratio_neg,fluid_name,number_of_test_points,sc_neg_soln,sc_pos_soln,Solvent_bead_SRD_box_density_cp_1):
    # Sc vs box size/lengthscale
    fig=plt.figure(figsize=(15,6))
    gs=GridSpec(nrows=1,ncols=1)
    fontsize=30

    #fig.suptitle(fluid_name+': $Sc\ vs$ $\\frac{\Delta x}{\\bar{\ell}}\\ $',size='large', wrap=True)

    ax1= fig.add_subplot(gs[0]) 
    for z in range(0,number_of_test_points):
        
        ax1.plot(mean_free_path_to_box_ratio_neg[z,:],sc_neg_soln[z,:],label='$M$={}'.format(int(Solvent_bead_SRD_box_density_cp_1[0,z])),marker='x')
        
        #ax1.legend(Solvent_bead_SRD_box_density_cp_1[0,z])
       # ax1.plot(mean_free_path_to_box_ratio_pos[z,:],sc_pos_soln[z,:],label='M+={}'.format(Solvent_bead_SRD_box_density_cp_1[0,z]))
        
        
        ax1.set_xscale('linear')
        ax1.set_yscale('log')
        ax1.set_xlabel('$Kn\  [-]$', rotation='horizontal',ha='right',fontsize=fontsize)
        ax1.set_ylabel( '$Sc\ [-]$', rotation='horizontal',ha='right',fontsize=fontsize)
        ax1.grid('on')
        ax1.legend(loc='right',bbox_to_anchor=(1.25, 0.5))
        
    plt.show()
        

sc_vs_mfp_to_collision_cell(mean_free_path_to_box_ratio_pos,mean_free_path_to_box_ratio_neg,fluid_name,number_of_test_points,sc_neg_soln,sc_pos_soln,Solvent_bead_SRD_box_density_cp_1)

# %%
from MPCD_constraints_on_solutions import MPCD_constraints 

no_timesteps=5000000
srd_ratio_tolerance=5
max_particle_count =1500000
min_particle_count=500
count_passed_constraints_neg=[]
count_passed_constraints_pos=[]
locations_of_non_nan_neg=()
locations_of_non_nan_pos=()
MPCD_constraints(no_timesteps,min_particle_count,sc_neg_soln,sc_pos_soln,srd_ratio_tolerance,max_particle_count,number_SRD_particles_wrt_pf_cp_mthd_1_pos,number_SRD_particles_wrt_pf_cp_mthd_1_neg,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,mean_free_path_pf_SRD_particles_cp_mthd_1_pos,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z],Solvent_bead_SRD_box_density_cp_1,tolerance,SRD_box_size_wrt_solid_beads[z],comparison_pos[z],comparison_neg[z])
# %%
#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Start of Ar Calculations #####
#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#%% Argon Calculations 
tolerance=0.01
number_of_test_points =50
Solvent_bead_SRD_box_density_cp_1 = np.array([(np.linspace(3.5,2000,number_of_test_points))])
number_of_M_cp_1=Solvent_bead_SRD_box_density_cp_1.shape[1]
number_boxes_var=128
min_number_boxes_for_particle_size=23
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)

#determine side length of simulaton box 
r_particle =50e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)
# define all inputs for Argon 
rho_s = 1426.9#621 #kg/m^3
r_particle =50e-6 #m 
#T_cel=34.5 #celsius, chosen from paper above 
Temp_visc_multiplier= 0.000099
T_K=86.5* Temp_visc_multiplier#+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST=0.00029800* Temp_visc_multiplier #Pa s


eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
# calculating stokes number in fluid conditions for solid particle tests
Stokes_number=0.0001
Gamma_dot= 4.5*Stokes_number*eta_s_NIST/ (rho_particle * r_particle**2)
fluid_name='Ar'

box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1.T

#%% dimensional SRD calculations 
from SRD_master_calc_dimensional import * 

SRD_dimensional_data= SRD_MASTER_calc_dimensional(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,T_K,k_b,rho_s,eta_s)
 
SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_dimensional_data[0]
SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_dimensional_data[1]
sc_pos_soln=SRD_dimensional_data[2]
sc_neg_soln=SRD_dimensional_data[3]
mean_free_path_pf_SRD_particles_cp_mthd_1_neg=SRD_dimensional_data[4]
mean_free_path_pf_SRD_particles_cp_mthd_1_pos=SRD_dimensional_data[5]
number_SRD_particles_wrt_pf_cp_mthd_1_neg=SRD_dimensional_data[6]
number_SRD_particles_wrt_pf_cp_mthd_1_pos=SRD_dimensional_data[7]
mass_fluid_particle_wrt_pf_cp_mthd_1=SRD_dimensional_data[8]
 
 
 
 
scaled_timestep=0.1
lengthscale=r_particle * 0.01
r_particle_scaled=r_particle/lengthscale

box_size_vec_nd= box_size_vec/lengthscale
mass_scale= 100000*rho_s* (lengthscale**3)
mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/mass_scale
energy_scale=k_b*T_K
timescale= np.sqrt((mass_scale*(lengthscale**2))/energy_scale)
scaled_temp=T_K/(energy_scale/k_b)
SRD_timestep_cp_1_based_on_sphere_neg_nd= SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check/timescale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=mean_free_path_pf_SRD_particles_cp_mthd_1_neg/lengthscale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_neg/box_size_vec
mean_free_path_pf_SRD_particles_cp_mthd_1_pos_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_pos/box_size_vec
SRD_MD_Ratio_neg= SRD_timestep_cp_1_based_on_sphere_neg_nd/scaled_timestep


#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Start of Water Calculations #####
#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#%%
fluid_name='Water'
tolerance=0.01
number_of_test_points =50
Solvent_bead_SRD_box_density_cp_1 = np.array([(np.linspace(3.5,2000,number_of_test_points))])
number_of_M_cp_1=Solvent_bead_SRD_box_density_cp_1.shape[1]
number_boxes_var=128
min_number_boxes_for_particle_size=23
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)

#determine side length of simulaton box 
r_particle =50e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)
rho_s = 1005##kg/m^3
r_particle =50e-6 #m 
#Temp_visc_multiplier=0.00000000099
Temp_visc_multiplier=0.0000001999#099 
T_K=300* Temp_visc_multiplier#+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST=0.00085253*Temp_visc_multiplier

#determine side length of simulaton box 
r_particle =50e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)
eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres

mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
# calculating stokes number in fluid conditions for solid particle tests
Stokes_number=0.0001
Gamma_dot= 4.5*Stokes_number*eta_s_NIST/ (rho_particle * r_particle**2)
box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1.T


# dimensional SRD calculations 
from SRD_master_calc_dimensional import * 

SRD_dimensional_data= SRD_MASTER_calc_dimensional(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,T_K,k_b,rho_s,eta_s)
 
SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_dimensional_data[0]
SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_dimensional_data[1]
sc_pos_soln=SRD_dimensional_data[2]
sc_neg_soln=SRD_dimensional_data[3]
mean_free_path_pf_SRD_particles_cp_mthd_1_neg=SRD_dimensional_data[4]
mean_free_path_pf_SRD_particles_cp_mthd_1_pos=SRD_dimensional_data[5]
number_SRD_particles_wrt_pf_cp_mthd_1_neg=SRD_dimensional_data[6]
number_SRD_particles_wrt_pf_cp_mthd_1_pos=SRD_dimensional_data[7]
mass_fluid_particle_wrt_pf_cp_mthd_1=SRD_dimensional_data[8]
 
 
 
 
scaled_timestep=0.01
lengthscale=r_particle * 0.1
r_particle_scaled=r_particle/lengthscale

box_size_vec_nd= box_size_vec/lengthscale
mass_scale= 100*rho_s* (lengthscale**3)
mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/mass_scale
energy_scale=k_b*T_K
timescale= np.sqrt((mass_scale*(lengthscale**2))/energy_scale)
scaled_temp=T_K/(energy_scale/k_b)
SRD_timestep_cp_1_based_on_sphere_neg_nd= SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check/timescale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=mean_free_path_pf_SRD_particles_cp_mthd_1_neg/lengthscale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_neg/box_size_vec
mean_free_path_pf_SRD_particles_cp_mthd_1_pos_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_pos/box_size_vec
SRD_MD_Ratio_neg= SRD_timestep_cp_1_based_on_sphere_neg_nd/scaled_timestep

#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Start of Cyclohex  Calculations #####
#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#%%
tolerance=0.01
number_of_test_points =50
Solvent_bead_SRD_box_density_cp_1 = np.array([(np.linspace(3.5,2000,number_of_test_points))])
number_of_M_cp_1=Solvent_bead_SRD_box_density_cp_1.shape[1]
number_boxes_var=128
min_number_boxes_for_particle_size=23
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)


rho_s = 764.95 #kg/m^3
r_particle =50e-6 #m 
T_cel=34.5 #celsius, chosen from paper above 
Temp_visc_multiplier=0.0000099#0.0000001999
T_K=T_cel+273.15 * Temp_visc_multiplier #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST=0.00076285 #Pa s 


eta_s=eta_s_NIST * Temp_visc_multiplier#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
r_particle =50e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)


Stokes_number=0.0001
Gamma_dot= 4.5*Stokes_number*eta_s_NIST/ (rho_particle * r_particle**2)
box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1.T
fluid_name='C6H12'


# dimensional SRD calculations 
from SRD_master_calc_dimensional import * 

SRD_dimensional_data= SRD_MASTER_calc_dimensional(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,T_K,k_b,rho_s,eta_s)
 
SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_dimensional_data[0]
SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_dimensional_data[1]
sc_pos_soln=SRD_dimensional_data[2]
sc_neg_soln=SRD_dimensional_data[3]
mean_free_path_pf_SRD_particles_cp_mthd_1_neg=SRD_dimensional_data[4]
mean_free_path_pf_SRD_particles_cp_mthd_1_pos=SRD_dimensional_data[5]
number_SRD_particles_wrt_pf_cp_mthd_1_neg=SRD_dimensional_data[6]
number_SRD_particles_wrt_pf_cp_mthd_1_pos=SRD_dimensional_data[7]
mass_fluid_particle_wrt_pf_cp_mthd_1=SRD_dimensional_data[8]
 
 
 
 
scaled_timestep=0.01
lengthscale=r_particle * 0.00001
r_particle_scaled=r_particle/lengthscale

box_size_vec_nd= box_size_vec/lengthscale
mass_scale= 100000*rho_s* (lengthscale**3)
mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/mass_scale
energy_scale=k_b*T_K
timescale= np.sqrt((mass_scale*(lengthscale**2))/energy_scale)
scaled_temp=T_K/(energy_scale/k_b)
SRD_timestep_cp_1_based_on_sphere_neg_nd= SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check/timescale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=mean_free_path_pf_SRD_particles_cp_mthd_1_neg/lengthscale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_neg/box_size_vec
mean_free_path_pf_SRD_particles_cp_mthd_1_pos_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_pos/box_size_vec
SRD_MD_Ratio_neg= SRD_timestep_cp_1_based_on_sphere_neg_nd/scaled_timestep

#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Start of hexane  Calculations #####
#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#%% 
tolerance=0.01
number_of_test_points =50
Solvent_bead_SRD_box_density_cp_1 = np.array([(np.linspace(3.5,2000,number_of_test_points))])
number_of_M_cp_1=Solvent_bead_SRD_box_density_cp_1.shape[1]
number_boxes_var=128
min_number_boxes_for_particle_size=23
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)

rho_s = 700#621 #kg/m^3
r_particle =50e-6 #m 
#T_cel=34.5 #celsius, chosen from paper above 
Temp_visc_multiplier=0.00099
T_K=311* Temp_visc_multiplier#+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST= 0.00046729 * Temp_visc_multiplier	 #Pa s 

eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
eta_s=eta_s_NIST * Temp_visc_multiplier#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
r_particle =50e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)
# calculating stokes number in fluid conditions for solid particle tests
Stokes_number=0.0001
Gamma_dot= 4.5*Stokes_number*eta_s_NIST/ (rho_particle * r_particle**2)

box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1.T
fluid_name='C6H14'

# dimensional SRD calculations 
from SRD_master_calc_dimensional import * 

SRD_dimensional_data= SRD_MASTER_calc_dimensional(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,T_K,k_b,rho_s,eta_s)
 
SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_dimensional_data[0]
SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_dimensional_data[1]
sc_pos_soln=SRD_dimensional_data[2]
sc_neg_soln=SRD_dimensional_data[3]
mean_free_path_pf_SRD_particles_cp_mthd_1_neg=SRD_dimensional_data[4]
mean_free_path_pf_SRD_particles_cp_mthd_1_pos=SRD_dimensional_data[5]
number_SRD_particles_wrt_pf_cp_mthd_1_neg=SRD_dimensional_data[6]
number_SRD_particles_wrt_pf_cp_mthd_1_pos=SRD_dimensional_data[7]
mass_fluid_particle_wrt_pf_cp_mthd_1=SRD_dimensional_data[8]
 
 
 
 
scaled_timestep=0.01
lengthscale=r_particle * 0.00001
r_particle_scaled=r_particle/lengthscale

box_size_vec_nd= box_size_vec/lengthscale
mass_scale= 100000*rho_s* (lengthscale**3)
mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/mass_scale
energy_scale=k_b*T_K
timescale= np.sqrt((mass_scale*(lengthscale**2))/energy_scale)
scaled_temp=T_K/(energy_scale/k_b)
SRD_timestep_cp_1_based_on_sphere_neg_nd= SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check/timescale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=mean_free_path_pf_SRD_particles_cp_mthd_1_neg/lengthscale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_neg/box_size_vec
mean_free_path_pf_SRD_particles_cp_mthd_1_pos_to_box_size= mean_free_path_pf_SRD_particles_cp_mthd_1_pos/box_size_vec
SRD_MD_Ratio_neg= SRD_timestep_cp_1_based_on_sphere_neg_nd/scaled_timestep

# %%
