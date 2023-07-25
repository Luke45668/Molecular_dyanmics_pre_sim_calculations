#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 25 2023

This script will do the SRD calculations in dimensional form, this will then be used to vary the temperature and viscosity scale. 



@author: lukedebono
"""
#%%
import os
import numpy as np
import sigfig
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


tolerance=0.01
Solvent_bead_SRD_box_density_cp_1 =10
# defining collision cell size range 

#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Start of N2 Calculations #####
################################################################
#%%
fluid_name='Nitrogen'
number_boxes_var=100 # 300
min_number_boxes_for_particle_size=6
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)
r_particle =25e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)


### physical parameters for Nitrogen
# fixed values for nitrogen 
rho_s = 847 #kg/m^3
number_of_test_points=10
#Temp_visc_multiplier=np.array([np.logspace(-6.5,0,number_of_test_points)]).T#0.0000099 #0.00005
Temp_visc_multiplier=np.array([np.geomspace(0.0000068129,0.00068129,number_of_test_points)]).T#0.0000099 #0.00005
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
nu_s_before_scale=eta_s_NIST/rho_s
eta_s_multiplier=1
eta_s=eta_s_NIST*Temp_visc_multiplier #*1000 to convert kg to g

nu_s = (eta_s/rho_s) 
box_size_vec = np.array([box_side_length/number_boxes_vec])

mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1
scaled_timestep=0.01
lengthscale=r_particle * 0.001
r_particle_scaled=r_particle/lengthscale

box_size_vec_nd= box_size_vec/lengthscale
mass_scale= 10000000*rho_s* (lengthscale**3)
#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# Start of Ar Calculations #####
#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#%%
# define all inputs for Argon 
fluid_name='Ar'
number_boxes_var=100 # 300
min_number_boxes_for_particle_size=6
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)
r_particle =25e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)

rho_s = 1426.9#621 #kg/m^3
r_particle =50e-6 #m 
#T_cel=34.5 #celsius, chosen from paper above 
number_of_test_points=10
Temp_visc_multiplier=np.array([np.logspace(-9,0,number_of_test_points)]).T
Temp_visc_multiplier=np.array([np.geomspace(0.0000083768,0.00083768,number_of_test_points)]).T
T_K=86.5 * Temp_visc_multiplier#+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST=0.00029800 #Pa s


eta_s=eta_s_NIST* Temp_visc_multiplier#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
# calculating stokes number in fluid conditions for solid particle tests
Stokes_number=0.0001
Gamma_dot= 4.5*Stokes_number*eta_s_NIST/ (rho_particle * r_particle**2)

scaled_timestep=0.1

box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1
#length_multiplier=np.repeat(np.array([np.logspace(-3,-1.5,number_of_lengthscales)]).T,number_boxes_var,axis=1)
length_multiplier= 10e-2 # chosen to be the midle of the reange 
lengthscale_parameter = length_multiplier*r_particle

#%% Water calculation 

fluid_name='Water'
number_boxes_var=100 # 300
min_number_boxes_for_particle_size=8
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)
r_particle =25e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)
scaled_timestep=0.01
tolerance=0.01
number_of_test_points =10

#determine side length of simulaton box 

rho_s = 1005##kg/m^3
r_particle =50e-6 #m 
Temp_visc_multiplier=np.array([np.logspace(-6,-7,number_of_test_points)]).T#099 
Temp_visc_multiplier=np.array([np.linspace(1e-6,6e-6,number_of_test_points)]).T#099 
T_K=300* Temp_visc_multiplier#+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST=0.00085253*Temp_visc_multiplier

 

eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres

mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
# calculating stokes number in fluid conditions for solid particle tests
Stokes_number=0.0001
Gamma_dot= 4.5*Stokes_number*eta_s_NIST/ (rho_particle * r_particle**2)
box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1

#box_size_vec_nd=box_side_length_scaled/number_boxes_vec
#SRD_box_size_wrt_solid_beads=box_size_vec_nd
#SRD_box_size_wrt_solid_beads_check=box_size_vec
#%% Cyclohexane 
fluid_name='C6H12'
number_boxes_var=100 # 300
min_number_boxes_for_particle_size=4
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)
r_particle =25e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)
rho_s = 764.95 #kg/m^3
r_particle =50e-6 #m 
T_cel=34.5 #celsius, chosen from paper above 
Temp_visc_multiplier=np.array([np.geomspace(0.00022,0.01,number_of_test_points)]).T
T_K=T_cel+273.15 * Temp_visc_multiplier #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST=0.00076285*Temp_visc_multiplier #Pa s 


eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1

length_multiplier=0.3162277#np.repeat(np.array([np.logspace(-1,0,number_of_lengthscales)]).T,number_boxes_var,axis=1)
#length_multiplier=np.repeat(np.array([np.logspace(-2.5,-1.5,number_of_lengthscales)]).T,number_boxes_var,axis=1)
lengthscale_parameter = length_multiplier*r_particle
box_side_length_scaled=(box_side_length/lengthscale_parameter)
box_size_to_lengthscale=box_size_vec/lengthscale_parameter
mass_multiplier=100
mass_scale = mass_multiplier* rho_s * (lengthscale_parameter**3)
r_particle_scaled = r_particle/lengthscale_parameter

box_size_vec = np.array([box_side_length/number_boxes_vec])

#%% Hexane 
fluid_name='C6H14'
number_boxes_var=100 # 300
min_number_boxes_for_particle_size=4
number_boxes_vec=np.linspace(min_number_boxes_for_particle_size,(min_number_boxes_for_particle_size-1)+number_boxes_var,number_boxes_var)
r_particle =25e-6
phi=0.005
N=2
Vol_box_at_specified_phi= N* (4/3)*np.pi*r_particle**3 /phi
box_side_length=np.cbrt(Vol_box_at_specified_phi)
rho_s = 700#621 #kg/m^3
r_particle =50e-6 #m 
#T_cel=34.5 #celsius, chosen from paper above 
Temp_visc_multiplier=np.array([np.logspace(-6,0,number_of_test_points)]).T
T_K=311* Temp_visc_multiplier#+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST= 0.00046729*Temp_visc_multiplier	 #Pa s 

eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
# calculating stokes number in fluid conditions for solid particle tests
Stokes_number=0.0001
Gamma_dot= 4.5*Stokes_number*eta_s_NIST/ (rho_particle * r_particle**2)
box_size_vec = np.array([box_side_length/number_boxes_vec])
mass_fluid_particle_wrt_pf_cp_mthd_1=(rho_s * (box_size_vec**3))/Solvent_bead_SRD_box_density_cp_1


#%% 





from SRD_master_calc_dimensional_not_changing_M import * 

SRD_dimensional_data= SRD_MASTER_calc_dimensional_no_M_variation(mass_fluid_particle_wrt_pf_cp_mthd_1,box_side_length,number_boxes_vec,tolerance,nu_s,Solvent_bead_SRD_box_density_cp_1 ,r_particle, box_size_vec ,T_K,k_b,rho_s,eta_s)
 
SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check=SRD_dimensional_data[0]
SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check=SRD_dimensional_data[1]
sc_pos_soln=SRD_dimensional_data[2]
sc_neg_soln=SRD_dimensional_data[3]
mean_free_path_pf_SRD_particles_cp_mthd_1_neg=SRD_dimensional_data[4]
mean_free_path_pf_SRD_particles_cp_mthd_1_pos=SRD_dimensional_data[5]
number_SRD_particles_wrt_pf_cp_mthd_1_neg=SRD_dimensional_data[6]
number_SRD_particles_wrt_pf_cp_mthd_1_pos=SRD_dimensional_data[7]
mass_fluid_particle_wrt_pf_cp_mthd_1=SRD_dimensional_data[8]

# need to change the length and mass scale for each fluid


mass_fluid_particle_wrt_pf_cp_mthd_1_nd=mass_fluid_particle_wrt_pf_cp_mthd_1/mass_scale
energy_scale=k_b*T_K
timescale= np.sqrt((mass_scale*(lengthscale**2))/energy_scale)
scaled_temp=T_K/(energy_scale/k_b)
SRD_timestep_cp_1_based_on_sphere_neg_nd= SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check/timescale
mean_free_path_pf_SRD_particles_cp_mthd_1_neg_nd=mean_free_path_pf_SRD_particles_cp_mthd_1_neg/lengthscale
box_size_vec_nd= box_size_vec/lengthscale
mean_free_path_to_box_ratio_neg= mean_free_path_pf_SRD_particles_cp_mthd_1_neg/box_size_vec
mean_free_path_to_box_ratio_pos= mean_free_path_pf_SRD_particles_cp_mthd_1_pos/box_size_vec
SRD_MD_Ratio_neg= SRD_timestep_cp_1_based_on_sphere_neg_nd/scaled_timestep

#%% plot Sc vs kn for various temp visc scalings
plt.rcParams.update({'font.size': 25})
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
fig=plt.figure(figsize=(15,6))
gs=GridSpec(nrows=1,ncols=1)
fontsize=20

#fig.suptitle(fluid_name+': $Sc\ vs$ $Kn\\ $',fontsize=fontsize, wrap=True)

ax1= fig.add_subplot(gs[0]) 
for z in range(0,number_of_test_points):
        
        ax1.plot(mean_free_path_to_box_ratio_neg[z,:],sc_neg_soln[z,:],label='$t_s$={}'.format((sigfig.round(Temp_visc_multiplier[z,0], sigfigs=2))),marker='x')
       
        #ax1.legend(Solvent_bead_SRD_box_density_cp_1[0,z])
       # ax1.plot(mean_free_path_to_box_ratio_pos[z,:],sc_pos_soln[z,:],label='M+={}'.format(Solvent_bead_SRD_box_density_cp_1[0,z]))
        
        
        ax1.set_xscale('linear')
        ax1.set_yscale('log')
        ax1.set_xlabel('$Kn\  [-]$', rotation='horizontal',ha='right')#,fontsize=fontsize)
        ax1.set_ylabel( '$Sc\ [-]$', rotation='horizontal',ha='right')#,fontsize=fontsize)
        ax1.grid('on')
        ax1.legend(loc='right',bbox_to_anchor=(1.3, 0.5))
        
plt.show()

#%% plot Kn vs collision cell size for various temp visc scalings
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
fig=plt.figure(figsize=(20,6))
gs=GridSpec(nrows=1,ncols=1)

fig.suptitle(fluid_name+': $Kn\ vs$ $\\frac{\Delta x}{\\bar{\ell}}\ $',size='large', wrap=True)

ax1= fig.add_subplot(gs[0]) 
for z in range(0,number_of_test_points):
        
        ax1.plot(box_size_vec_nd[0,:],mean_free_path_to_box_ratio_neg[z,:],label='$t_s$={}'.format(Temp_visc_multiplier[z,0]),marker='x')
        
        #ax1.legend(Solvent_bead_SRD_box_density_cp_1[0,z])
       # ax1.plot(mean_free_path_to_box_ratio_pos[z,:],sc_pos_soln[z,:],label='M+={}'.format(Solvent_bead_SRD_box_density_cp_1[0,z]))
        
        
        ax1.set_xscale('linear')
        ax1.set_yscale('log')
        ax1.set_xlabel('$\\frac{\Delta x}{\\bar{\ell}}\  [-]$', rotation='horizontal',ha='right',size='large')
        ax1.set_ylabel( '$Kn\ [-]$', rotation='horizontal',ha='right',size='large')
        ax1.grid('on')
        ax1.legend(loc='right',bbox_to_anchor=(1.1, 0.5))
        
plt.show()
 

#%% plot Sc vs collision cell size for various temp visc scalings
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
fig=plt.figure(figsize=(20,6))
gs=GridSpec(nrows=1,ncols=1)

fig.suptitle(fluid_name+': $Sc\ vs$ $\\frac{\Delta x}{\\bar{\ell}}\ $',size='large', wrap=True)

ax1= fig.add_subplot(gs[0]) 
for z in range(0,number_of_test_points):
        
        ax1.plot(box_size_vec_nd[0,:],sc_neg_soln[z,:],label='$t_s$={}'.format((sigfig.round(Temp_visc_multiplier[z,0], sigfigs=2))),marker='x')
        
        #ax1.legend(Solvent_bead_SRD_box_density_cp_1[0,z])
       # ax1.plot(mean_free_path_to_box_ratio_pos[z,:],sc_pos_soln[z,:],label='M+={}'.format(Solvent_bead_SRD_box_density_cp_1[0,z]))
        
        
        ax1.set_xscale('linear')
        ax1.set_yscale('log')
        ax1.set_xlabel('$\\frac{\Delta x}{\\bar{\ell}}\  [-]$', rotation='horizontal',ha='right',size='large')
        ax1.set_ylabel( '$Sc\ [-]$', rotation='horizontal',ha='right',size='large')
        ax1.grid('on')
        ax1.legend(loc='right',bbox_to_anchor=(1.1, 0.5))
        
plt.show()
 
#%% plot Delta_tsrd vs collision cell size for various temp visc scalings
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
fig=plt.figure(figsize=(15,6))
gs=GridSpec(nrows=1,ncols=1)
fontsize=20

#fig.suptitle(fluid_name+': $\\frac{\Delta t_{SRD}}{\\tau}\ vs$ $\\frac{\Delta x}{\ell}\ $',fontsize=fontsize, wrap=True)

ax1= fig.add_subplot(gs[0]) 
for z in range(0,number_of_test_points):
        
        ax1.plot(box_size_vec_nd[0,:],SRD_timestep_cp_1_based_on_sphere_neg_nd[z,:],label='$t_s$={}'.format((sigfig.round(Temp_visc_multiplier[z,0], sigfigs=2))),marker='x')
        
        #ax1.legend(Solvent_bead_SRD_box_density_cp_1[0,z])
       # ax1.plot(mean_free_path_to_box_ratio_pos[z,:],sc_pos_soln[z,:],label='M+={}'.format(Solvent_bead_SRD_box_density_cp_1[0,z]))
        
        
        ax1.set_xscale('linear')
        ax1.set_yscale('log')
        ax1.set_xlabel('$\\frac{\Delta x}{\ell}\  [\ell]$', rotation='horizontal',ha='right',fontsize=fontsize)
        ax1.set_ylabel('$\\frac{\Delta t_{SRD}}{\\tau}\ [\\tau] $', rotation='horizontal',ha='right',fontsize=fontsize)
        ax1.grid('on')
        ax1.legend(loc='right',bbox_to_anchor=(1.1, 0.5))
        
plt.show()

#%% plot SRD MD vs collision cell size for various temp visc scalings
plt.rcParams.update({'font.size': 25})
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
fig=plt.figure(figsize=(15,6))
fontsize=20
gs=GridSpec(nrows=1,ncols=1)

#fig.suptitle(fluid_name+': $\\frac{\Delta t_{SRD}}{\Delta t_{MD}}\ vs$ $\\frac{\Delta x}{\\bar{\ell}}\ $',size='x-large', wrap=True)

ax1= fig.add_subplot(gs[0]) 
for z in range(0,number_of_test_points):
        
        ax1.plot(box_size_vec_nd[0,:],SRD_MD_Ratio_neg[z,:],label='$t_s$={}'.format((sigfig.round(Temp_visc_multiplier[z,0], sigfigs=2))),marker='x')
        
        #ax1.legend(Solvent_bead_SRD_box_density_cp_1[0,z])
       # ax1.plot(mean_free_path_to_box_ratio_pos[z,:],sc_pos_soln[z,:],label='M+={}'.format(Solvent_bead_SRD_box_density_cp_1[0,z]))
        
        
        ax1.set_xscale('linear')
        ax1.set_yscale('log')
        ax1.set_xlabel('$\\frac{\Delta x}{\\bar{\ell}}\  [\\tau]$', rotation='horizontal',ha='right')#fontsize=fontsize)
        ax1.set_ylabel('$\\frac{\Delta t_{SRD}}{\Delta t_{MD}}\ [-]$', rotation='horizontal',ha='right')#fontsize=fontsize)
        ax1.grid('on')
        ax1.legend(loc='best',bbox_to_anchor=(1.3, 1.05))
        
plt.show()

# %%
