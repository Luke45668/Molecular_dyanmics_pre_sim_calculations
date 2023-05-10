#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 09:53:07 2022

@author: lukedebono
"""

#import math as m
import os
import numpy as np
import numba
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




#%% Required constants

rho_s = 764.95 #kg/m^3
r_particle =50e-6 #m 
T_cel=34.5 #celsius, chosen from paper above 
T_K=T_cel+273.15 #Kelvin
k_b= 1.380649e-23 #boltzmann in J K^-1
eta_s_NIST=0.00076285 #Pa s or kg m^-1 s^-1 this is cyclohexane as our theta solvent at 34.5 degrees and 1atm
# Cyclohexane was chosen based on Jones, G., &#38; Caroline, D. (1978). 
#Intramolecular motion of polystyrene in a theta solvent.
# https://doi.org/10.1016/0009-2614(78)80336-4</div>
eta_s=eta_s_NIST#*1000 #*1000 to convert kg to g
nu_s = eta_s/rho_s
rho_particle = 1200 #kg m^-3 PMMA spheres
Stokes_number=0.0001
eta_polymer_fluid= np.array([[0.07682105, 0.04892482, 0.03435542]])
Gamma_dot= 4.5*Stokes_number*eta_polymer_fluid/ (rho_particle * r_particle**2)
number_of_test_points =30000
Solvent_bead_SRD_box_density_cp_1 = np.array([np.linspace(3.5,1000,number_of_test_points)])
number_of_M_cp_1=Solvent_bead_SRD_box_density_cp_1.shape[1]

ms=5 #plot marker size 
scaled_timestep = 0.01 
peclet_number_dimensional =( 6 * np.pi * (r_particle**3) * Gamma_dot * eta_s)/k_b * T_K
print(peclet_number_dimensional)
#%% Non-dimensionalising inputs directly from lammps
# fundamental scalings 
mass_solid_particle= rho_particle * (4/3)*np.pi*(r_particle**3)
SRD_mass_scale_parameter = mass_solid_particle
lengthscale_parameter = r_particle 
r_particle_scaled = r_particle/lengthscale_parameter
timescale_parameter=(1/Gamma_dot[0,1])*0.01
scaled_temp =5e-05
#%% determiningfundamental scalings
import units_lj_scalings
energy_parameter=units_lj_scalings.units_lj_scalings(SRD_mass_scale_parameter,lengthscale_parameter,timescale_parameter,k_b,rho_s,eta_s)[0]
timescale_parameter=units_lj_scalings.units_lj_scalings(SRD_mass_scale_parameter,lengthscale_parameter,timescale_parameter,k_b,rho_s,eta_s)[1]
temperature_parameter=units_lj_scalings.units_lj_scalings(SRD_mass_scale_parameter,lengthscale_parameter,timescale_parameter,k_b,rho_s,eta_s)[2]
scaled_dynamic_viscosity=units_lj_scalings.units_lj_scalings(SRD_mass_scale_parameter,lengthscale_parameter,timescale_parameter,k_b,rho_s,eta_s)[3]
scaled_nu_s=units_lj_scalings.units_lj_scalings(SRD_mass_scale_parameter,lengthscale_parameter,timescale_parameter,k_b,rho_s,eta_s)[4]
scaled_rho_s=units_lj_scalings.units_lj_scalings(SRD_mass_scale_parameter,lengthscale_parameter,timescale_parameter,k_b,rho_s,eta_s)[5]
#%%  possible SRD bin sizes
 
new_simbox_side_length = r_particle_scaled*16 
box_side_length_scaled = new_simbox_side_length
box_size_vec =np.array([[0.25,0.4,0.5,0.8,1,2,4]]) 
SRD_box_size_wrt_solid_beads = r_particle_scaled *box_size_vec 
SRD_box_size_wrt_solid_beads_check = r_particle *box_size_vec  

# gamma_1 = (1 - ((1-np.exp(-Solvent_bead_SRD_box_density_cp_1))/Solvent_bead_SRD_box_density_cp_1)).T
# gamma_2 = (Solvent_bead_SRD_box_density_cp_1* (Solvent_bead_SRD_box_density_cp_1+2)/(Solvent_bead_SRD_box_density_cp_1-1)).T  
import SRD_a_minus_gamma_funcs

gamma_1=SRD_a_minus_gamma_funcs.SRD_a_minus_gamma_funcs(Solvent_bead_SRD_box_density_cp_1)[0]
gamma_2=SRD_a_minus_gamma_funcs.SRD_a_minus_gamma_funcs(Solvent_bead_SRD_box_density_cp_1)[1]


#%% consraints for box size and derived constraint

import box_size_dim_and_integer_bin_constraint
box_size_dim_and_integer_bin_constraint.box_size_dim_and_integer_SRD_bin_count_constraint_func(k_b,T_K,gamma_1,rho_s,nu_s,box_side_length_scaled,SRD_box_size_wrt_solid_beads,SRD_box_size_wrt_solid_beads_check)

#%% Calculatiung mass of MPCD particles

volume_solid_particle = (4/3)*np.pi*(r_particle_scaled**3)
volume_solvent_without_solid = new_simbox_side_length**(3)

Number_of_SRD_boxes_in_sim_box_wrt_pf = (new_simbox_side_length**3)/(SRD_box_size_wrt_solid_beads**3)
number_SRD_particles_wrt_pf_cp_mthd_1 = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf * Solvent_bead_SRD_box_density_cp_1.T)
mass_fluid_particle_wrt_pf_cp_mthd_1= (volume_solvent_without_solid*scaled_rho_s)/number_SRD_particles_wrt_pf_cp_mthd_1

#%% calculating SRD timestep non-dimensional


import srd_timestep_non_dim
SRD_timestep_cp_1_based_on_sphere_pf_neg_nd = srd_timestep_non_dim.srd_timestep_non_dim(scaled_nu_s, scaled_rho_s, energy_parameter, timescale_parameter, SRD_mass_scale_parameter, lengthscale_parameter, scaled_temp, gamma_1, gamma_2, SRD_box_size_wrt_solid_beads)[0]
SRD_timestep_cp_1_based_on_sphere_pf_pos_nd=srd_timestep_non_dim.srd_timestep_non_dim(scaled_nu_s, scaled_rho_s, energy_parameter, timescale_parameter, SRD_mass_scale_parameter, lengthscale_parameter, scaled_temp, gamma_1, gamma_2, SRD_box_size_wrt_solid_beads)[1]

#print(SRD_timestep_cp_1_based_on_sphere_pf_neg_nd)

#%% testing the dimensional calculation against the N_d calc

import srd_timestep_dim
SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check =  srd_timestep_dim.srd_timestep_dim(nu_s,k_b,T_K,gamma_1,gamma_2,rho_s,SRD_box_size_wrt_solid_beads_check,timescale_parameter)[0]
SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check =  srd_timestep_dim.srd_timestep_dim(nu_s,k_b,T_K,gamma_1,gamma_2,rho_s,SRD_box_size_wrt_solid_beads_check,timescale_parameter)[1]
atol=0.01
rtol=0.1
arr=SRD_timestep_cp_1_based_on_sphere_pf_neg_nd 
arr1=SRD_timestep_cp_1_based_on_sphere_neg_dimensional_check
print("Negative n-d solution matches d solution is T/F?")
print((np.allclose(arr,arr1,rtol=rtol, equal_nan=True)))
arr=SRD_timestep_cp_1_based_on_sphere_pf_pos_nd 
arr1=SRD_timestep_cp_1_based_on_sphere_pos_dimensional_check
print("Positve n-d solution matches d solution is T/F?")
print((np.allclose(arr,arr1,rtol=rtol, equal_nan=True)))
       
     
#%% mean free path and timestep ratio 
mean_free_path_pf_SRD_particles_cp_mthd_1_pos =SRD_timestep_cp_1_based_on_sphere_pf_pos_nd* np.sqrt((scaled_temp)/(mass_fluid_particle_wrt_pf_cp_mthd_1))
mean_free_path_pf_SRD_particles_cp_mthd_1_neg = SRD_timestep_cp_1_based_on_sphere_pf_neg_nd* np.sqrt((scaled_temp)/(mass_fluid_particle_wrt_pf_cp_mthd_1))

number_SRD_particles_wrt_pf_cp_mthd_1_pos = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf * Solvent_bead_SRD_box_density_cp_1.T)
number_SRD_particles_wrt_pf_cp_mthd_1_neg = np.ceil(Number_of_SRD_boxes_in_sim_box_wrt_pf * Solvent_bead_SRD_box_density_cp_1.T)

Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos = SRD_timestep_cp_1_based_on_sphere_pf_pos_nd/scaled_timestep
Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg = SRD_timestep_cp_1_based_on_sphere_pf_neg_nd/scaled_timestep
#%% removing timestep ratios which arent close to integer values

check=np.zeros((Solvent_bead_SRD_box_density_cp_1.size,box_size_vec.size))
check= check.astype('float64')
comparison_pos = np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos)-Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos
comparison_neg = np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)-Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg

#%% checking for integer srd/md ratio and knudsen number continuum constraint 
tolerance=0.01
max_particle_count=100000
import MPCD_constraints_on_solutions

MPCD_constraints_on_solutions.MPCD_constraints_on_solutions(max_particle_count,number_SRD_particles_wrt_pf_cp_mthd_1_pos,number_SRD_particles_wrt_pf_cp_mthd_1_neg,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,mean_free_path_pf_SRD_particles_cp_mthd_1_pos,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,Solvent_bead_SRD_box_density_cp_1,tolerance,SRD_box_size_wrt_solid_beads,comparison_pos,comparison_neg) 

             
#%% creating a matrix of inputs
count_passed_constraints_neg = np.count_nonzero(~np.isnan(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)) 
# this counts the non-nan values of the array by inverting the true false routine of .isnan with a ~ so now false are 1s and trues are 0 
count_passed_constraints_pos = np.count_nonzero(~np.isnan(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos)) 
locations_of_non_nan_neg= np.argwhere(~np.isnan(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg))
locations_of_non_nan_pos= np.argwhere(~np.isnan(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos))

#%%

#### plotting dtSRD/dtMD vs alpha 
# Use the plot from this box
# using grid specs to plot 
# organising results 
x_data_cp_1 = Solvent_bead_SRD_box_density_cp_1.flatten()
#print(x_data_cp_1)
y_data_cp_1 = np.repeat(SRD_box_size_wrt_solid_beads,number_of_M_cp_1,axis=0)
#print(y_data_cp_1)
z_data_cp_1_pos = mean_free_path_pf_SRD_particles_cp_mthd_1_pos
z_data_cp_1_neg = mean_free_path_pf_SRD_particles_cp_mthd_1_neg

fig=plt.figure(figsize=(10,6))
gs=GridSpec(nrows=1,ncols=1)
# ax1= fig.add_subplot(1,1,1,adjustable='box', aspect=1)
# ax2= fig.add_subplot(1,2,2,adjustable='box', aspect=1)
#fig, (ax1,ax2)=plt.subplots(1,2)
fig.suptitle('SRD Solid-Liquid Coupling Method 1(not inc. Solid sphere): Mean free path vs Collision cell size $\delta t_{MD}='+str(scaled_timestep)+'$',size='large',ha='center')

ax1= fig.add_subplot(gs[0]) #[lowerCorner_x, lowerCorner_y, width, height]
#ax2= fig.add_subplot(gs[1])     #[lowerCorner_x, lowerCorner_y, width, height]

# if you need to reset stle as plt.style is not responding use the line below
plt.rcParams.update(plt.rcParamsDefault) 
#seaborn-notebook
#plt.style.use('seaborn')
#plt.subplots_adjust(hspace=5)
for i in range (0,number_of_M_cp_1):
    #for j in range (0,4): 
     ax1.plot( y_data_cp_1[i,:], z_data_cp_1_pos[i,:],ms=ms,label='M+={}'.format(x_data_cp_1[i]),linestyle='--', marker="x")
     ax1.plot( y_data_cp_1[i,:], z_data_cp_1_neg[i,:],ms=ms,label='M-={}'.format(x_data_cp_1[i]),linestyle='--', marker="+")
     ax1.legend(x_data_cp_1)
     ax1.set_xscale('linear')
     ax1.set_yscale('log')
     #ax1.set_title('SRD Coupling Method 1', size='x-large')
     ax1.set_xlabel('$\Delta x^{*} [-]$',size='large')
     ax1.set_ylabel( '$\lambda^{*} [-]$', rotation='horizontal',ha='right',size='large')
     ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left',borderaxespad=1)
    
     ax1.grid()
    # ax1.imshow( aspect='auto')
#on y asix
     
# =============================================================================
# for i in range (0,number_of_M_cp_2):
#     #for j in range (0,4): 
#      ax2.plot( y_data_cp_2[i,:], z_data_cp_2[i,:],label='M={}'.format(x_data_cp_2[i]),linestyle='--', marker='o')
#      ax2.legend(x_data_cp_1)
#    
#      ax2.set_xscale('log')
#      ax2.grid()
#      ax2.set_yscale('log')
#      ax2.set_title('SRD Coupling Method 2', size='x-large') 
#      ax2.set_xlabel('$\lambda/a$',size='x-large')
#      ax2.set_ylabel( '$\delta t_{SRD}/\delta t_{MD}$', rotation='horizontal',ha='right',size='x-large')
#      ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=1)
#      #ax2.set_aspect(aspect='equal',adjustable='box')
# 
#      
# =============================================================================
fig.tight_layout()
# Most likely, the problem is that you're using a relative file path to open the file, but the current working directory isn't set to what you think it is.

# It's a common misconception that relative paths are relative to the location of the python script, but this is untrue. Relative file paths are always relative to the current working directory, and the current working directory doesn't have to be the location of your python script.

# You have three options:

# Use an absolute path to open the file.
# Generate the path to the file relative to your python script.
# Change the current working directory before opening the file
#path= os.chdir("/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Plots ")
#plt.savefig("mean_free_path_vs_collision_cell_size_not_inc_sphere_timestep_no_pos_data_M_2.5_10_dt_"+str(scaled_timestep)+".png")
plt.show()
#os.chdir("/Volumes/Backup Plus/PhD_/Rouse Model simulations")
#os.getcwd()


# formatting keeps changing between runs and i've not changed the code at all.




#%% Simulation Run section cp 1
# for imac
#os.chdir('/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac') 
# for laptop
os.chdir('/Users/lukedebono/Documents/LAMMPS_projects_mac_book/OneDrive_1_24-02-2023/LAMMPS python run and analysis scripts/Analysis codes') 
import velP2numpy 
import log2numpy
l=10 # number of repeats
m=0 # alpha value selection 
k=0
#initialising the vars 
realisation_index_ =np.linspace(0, l,l) # [1,2,3]
stress_ave = np.zeros((l,15))
apparent_visc =np.zeros((l,15))
SRD_temp = np.zeros((l,15))
lamda=mean_free_path_pf_SRD_particles_cp_mthd_1_neg#np.zeros((l,20))
equilibration_timesteps= 1000 # number of steps to do equilibration with 
VP_ave_freq =1000
chunk = 20
#swap_rate=np.array([range(3,1202,3)])
swap_rate = np.array([3,7,15,30,60,300,600,900,1200]) # values chosen from original mp paper
swap_number = np.array([1,10,100,1000])
simulation_file="Simulation_run_folder"
number_of_realisations=5 
number_of_M_values = Solvent_bead_SRD_box_density_cp_1.size
shear_rate_average=np.zeros((swap_rate.size,number_of_M_values,number_of_realisations))# change 1 if you use various swap rates
no_timesteps=500000
dump_freq=1000 # if you change the timestep rememebr to chaneg this 
thermo_freq = 10000
path_2_simulation_run ='/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder'
#VP_data_upper=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size,9,int(no_timesteps/thermo_freq),5))
#VP_data_lower=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size,9,int(no_timesteps/thermo_freq),5))
#%% simulation run code for neg solns on home computer
import run_mpi_simulation_neg_soln
run_mpi_simulation_neg_soln.run_mpi_simulation_neg_soln(equilibration_timesteps,VP_ave_freq,realisation_index_,path_2_simulation_run,scaled_temp,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,lamda,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,scaled_timestep,box_side_length_scaled,mass_fluid_particle_wrt_pf_cp_mthd_1,number_SRD_particles_wrt_pf_cp_mthd_1,swap_number,swap_rate,locations_of_non_nan_neg,number_of_realisations)

# Path_2_VP="/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/"+simulation_file
# VP_output_col_count = 4 
# realisation_name="vel.profile_no_tstat__no_rescale_"+str(rand_int)+"_"+str(realisation_index_[j])+"_"+str(no_SRD)+"_"+str(box_size)+"_"+str(timestep_input)+"_"+str(SRD_MD_ratio)+"_"+str(dump_freq)+"_"+str(thermo_freq)+"_"+str(no_timesteps)+"_T_"+str(temp_)+"_lbda_"+str(lamda[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])+"_SR_"+str(swap_rate[m])+"_SN_"+str(swap_number[k])  
# VP_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[0]
# VP_z_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[1]
# #          # split the profile into mirror images  
# VP_data_upper[z,m,:,:,j] = VP_data[1:10,:]
# VP_z_data_upper = VP_z_data[1:10].astype('float64')* box_side_length_scaled
# VP_data_lower[z,m,:,:,j] = VP_data[11:,:]
# VP_z_data_lower = VP_z_data[11:].astype('float64') * box_side_length_scaled
#%% writing input data to a parameter file 
"""
Format will be the order the data points go into lammps, this can then be made into a general function that writes the file in the order you want 

Probably advisable to put key parameters and creation meta data in the name 

"""
#os.mkdir('Parameter_arrays_file')
os.chdir('/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Parameter_arrays_file')
array_entry = np.linspace(1,150,150)
name_lammps_script='no_wall_pure_SRD_sim_var_inputs_td_var_no_tstat_no_rescale_no_dump.file'
lammps_input_parameters_list = ' Array_entry swap_rate swap_number VP_ave_freq equilibration_timesteps chunk  grid_size realisation_index temp_ lambda rand_int rand_int_1 rand_int_2 no_SRD mass_SRD box_size timestep_input SRD_MD_ratio dump_freq thermo_freq no_timesteps'

for z in range(0,locations_of_non_nan_neg.shape[0]): 
    with open('input_parameters_solution_no_'+str(z)+'_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_no_swap_freqs_'+str(swap_rate.size)+'_no_test_points_'+str(number_of_test_points)+'.txt','w') as f:
        f.write('Parameter file for LAMMPS script: '+name_lammps_script+' \n Param order sequence: \n '+lammps_input_parameters_list+'\n')
  
       # for i in range (0,7): 
           # i=1
        for j in range(0,3): #or now just use one realisation 
           for m in range(0,swap_rate.size):#range(0,1):       #range(0,swap_rate.size):
                    for k in range(0,swap_number.size):
                        #or l in range(0,array_entry.size):
                        
                            no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])) 
                            mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])
                            box_size = str(box_side_length_scaled)
                            timestep_input= str(scaled_timestep)
                            chunk = 20 # number of chunks to use for VP averaging
                            SRD_MD_ratio=str(int(np.round( Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])))
                            lamda[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]] = mean_free_path_pf_SRD_particles_cp_mthd_1_neg[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]] 
                            grid_size = str(SRD_box_size_wrt_solid_beads[0,locations_of_non_nan_neg[z,1]])
                            dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                            thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                            no_timesteps = str(no_timesteps)
                            no_timesteps_ =float(no_timesteps)
                            temp_=str(scaled_temp)
                            rand_int =str(np.random.randint(0, 1000000))
                            rand_int_1 =str( np.random.randint(0, 1000000))
                            rand_int_2 =str(np.random.randint(0, 1000000))
                            f.write( str(int(array_entry[z]))+' '+str(swap_rate[m])+' '+ str(swap_number[k])+' '+ str(VP_ave_freq)+' '+ str(equilibration_timesteps)+' '+ str(chunk)+' '+ str(grid_size)+' '+ str(realisation_index_[j])+' '+ str(temp_)+' '+ str(lamda[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])+' '+ rand_int+' '+ rand_int_1+' '+ rand_int_2 +' '+no_SRD+' '+ mass_SRD+' '+ box_size+' '+ timestep_input+' '+ SRD_MD_ratio+' '+ dump_freq+' '+ thermo_freq+' '+ no_timesteps+'\n' )
    # basically need to put an f.write with all the vars in the right order in a for loop to write the whole sequence , also 
# need to produce the array entry list 
#%% Fixed overwrite variables
"""
This section is option 2 for production runs, write all the realisations for each solution to a single bash script. 


"""
os.chdir('/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac')
META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
j_=3
Path_2_generic='/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD/'
specific_email = 'luke.debono.21@ucl.ac.uk'
simulation_batch_folder= '/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD/simulation_batch_scripts_validations_realisation_'+str(j_)+'_'+'_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_no_swap_freqs_'+str(swap_rate.size)+'_no_test_points_'+str(number_of_test_points)+'_'+META_DATA 
os.mkdir(simulation_batch_folder)
#os.chdir('/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD/simulation_run_val_'+'_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_no_swap_freqs_'+str(swap_rate.size)+'_no_test_points_'+str(number_of_test_points)+META_DATA )
wall_time='4:00:00'
ram_requirement='4G'# per task 
tempdir_req='50G'
num_task_req=''
np=4
# repeat number
total_no_realisations_per_solution=9 
np_req=str(total_no_realisations_per_solution*np)
wd_path='/home/ucahlrl/Scratch/output/' #simulation_batch_validations_'+'_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_no_swap_freqs_'+str(swap_rate.size)+'_no_test_points_'+str(number_of_test_points)+'_'+META_DATA 
extra_code='module load mpi/intel/2018/update3/intel\n'
data_transfer_instructions=''



# for loop to create run line 

#%% creating scripts for qsub 
import py2bash_launch_overwriter
import numpy as np

for z in range(0,locations_of_non_nan_neg.shape[0]): 
      #run_code=''
      #simulation_run_name='cyclohexane_val_param_sweep_solution_'+str(z)+'_test_run_'
      for j in range(j_-1,j_): #or now just use one realisation 
       
           for k in range(0,swap_number.size):
                   simulation_run_name='cyclohexane_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_number_'+str(swap_number[k])+'_test_run_'
                   run_code=''   
                   for m in range(0,swap_rate.size):#range(0,1):  
                        #or l in range(0,array_entry.size):
                        
                            no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])) 
                            print(no_SRD)
                            mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])
                            box_size = str(box_side_length_scaled)
                            timestep_input= str(scaled_timestep)
                            chunk = 20 # number of chunks to use for VP averaging
                            SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])
                            SRD_MD_ratio=str(int(SRD_MD_ratio))
                            lamda[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]] = mean_free_path_pf_SRD_particles_cp_mthd_1_neg[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]] 
                            grid_size = str(SRD_box_size_wrt_solid_beads[0,locations_of_non_nan_neg[z,1]])
                            dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                            thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                            no_timesteps = str(no_timesteps)
                            no_timesteps_ =float(no_timesteps)
                            temp_=str(scaled_temp)
                            rand_int =str(np.random.randint(0, 1000000))
                            rand_int_1 =str( np.random.randint(0, 1000000))
                            rand_int_2 =str(np.random.randint(0, 1000000))


                            run_code_individual ='mpirun -np 4 /home/ucahlrl/simulation_run_folder/lammps-23Jun2022/src/lmp_mpi -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]])+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+'  -in /home/ucahlrl/simulation_run_folder/no_wall_pure_SRD_sim_var_inputs_td_var_no_tstat_no_rescale_no_dump.file & \n'
                            run_code=run_code +run_code_individual

                   run_code = run_code[:-3]

                   py2bash_launch_overwriter.py2bash_launch_overwriter(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)

#%% save the vp arrays                     
# os.chdir("/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder")
# np.save("VP_data_upper_T"+str(scaled_temp)+"_NoTP_"+str(number_of_test_points)+"no_timesteps_"+str(no_timesteps),VP_data_upper)
# np.save("VP_data_lower_T"+str(scaled_temp)+"_NoTP_"+str(number_of_test_points)+"no_timesteps_"+str(no_timesteps),VP_data_lower)
#%% reading in old data
import glob 
VP_output_col_count = 4 

count=0
realisation_name = [""]     
#for name in glob.glob('/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder/VP_data using pearson coeff/vel.profile_no_tstat__no_rescale_*'):
# make sure you are in the correct location before you read in the data
for name in glob.glob('vel.profile_no_tstat__no_rescale_*'):
    count=count+1    
    realisation_name.append(name)
    print(name)
# getting rid of the first empty cell    
realisation_name =realisation_name[1:]     

#%% organising data after reading in 
VP_output_col_count = 4 
count_1=0
simulation_file="Simulation_run_folder/VP_data using pearson coeff"
Path_2_VP="/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/"+simulation_file

# find the realisation names that match 
for z in range(0,locations_of_non_nan_neg.shape[0]): 
   # for i in range (0,7): 
       # i=1
        for j in range(0,1):# for now just use one realisation 
            for m in range(0,swap_rate.size):#range(0,1):       #range(0,swap_rate.size):
                for k in range(0,1):#,swap_number.size):
                    # find the realisation names that match 
                    no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]]))
                    glob_search='vel.profile_no_tstat__no_rescale_*'+no_SRD+"*_SR_"+str(swap_rate[m])+"_SN_"+str(swap_number[k])
                     
                    for name in glob.glob(glob_search):
                        realisation_name=name 
                        
                        VP_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[0]
                        VP_z_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[1]
                    #          # split the profile into mirror images  
                        VP_data_upper[z,m,:,:,j] = VP_data[1:10,:]
                        VP_z_data_upper = VP_z_data[1:10].astype('float64')* box_side_length_scaled
                        VP_data_lower[z,m,:,:,j] = VP_data[11:,:]
                        VP_z_data_lower = VP_z_data[11:].astype('float64') * box_side_length_scaled

#%% Regression line testing code 
# importing the data 
VP_data_upper=VP_data_upper_T5e05_NoTP_10000no_timesteps_200000
VP_data_lower=VP_data_lower_T5e05_NoTP_10000no_timesteps_200000
# organising the data

VP_data_lower_realisation_averaged = np.mean(VP_data_lower,axis=4) 
# 4 indicates to take the mean of the 
VP_data_upper_realisation_averaged = np.mean(VP_data_upper,axis=4) 


#%% plotting pearson coeff vs 
# upper VP

x= np.array(VP_data_upper_realisation_averaged)
y= np.repeat(np.array([VP_z_data_upper]).T,VP_data_upper.shape[3],1)
pearson_coeff_upper= np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size,VP_data_upper_realisation_averaged.shape[3]))
shear_rate_upper= np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size,VP_data_upper_realisation_averaged.shape[3]))
shear_rate_grad_upper=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size,VP_data_upper_realisation_averaged.shape[3]))
#pearson_mean_upper=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size))    
#standard_deviation_upper=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size))                  
perfect_linearity_comparison =np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size,VP_data_upper_realisation_averaged.shape[3]))
t_s=np.array([[np.linspace(0,VP_data_upper.shape[3],int(float(no_timesteps)/float(thermo_freq)))]])
t_s=np.repeat(t_s, locations_of_non_nan_neg.shape[0],axis=0)
t_s=np.repeat(t_s, swap_rate.size,axis=1)     
swap_rate_for_plotting_ = np.array([[swap_rate]])    
swap_rate_for_plotting_=np.repeat(swap_rate_for_plotting_,locations_of_non_nan_neg.shape[0],axis=0)
swap_rate_for_plotting_=np.repeat(swap_rate_for_plotting_,int(float(no_timesteps)/float(thermo_freq)),axis=1)

#%%                              
for z in range(0,locations_of_non_nan_neg.shape[0]): 
    for m in range(0,swap_rate.size):
        for i in range(0,VP_data_upper.shape[3]):
                
                pearson_coeff_upper[z,m,i] =scipy.stats.pearsonr(x[z,m,:,i], y[:,i])[0]
                shear_rate_upper[z,m,i]= scipy.stats.linregress(y[:,i],x[z,m,:,i] ).slope
               
                
                #pearson_mean_upper[z,m]= np.mean(pearson_coeff_upper[z,m,:])
                perfect_linearity_comparison[z,m,i]=  1-pearson_coeff_upper[z,m,i]
                # standard_deviation_upper [z,m]=np.std(pearson_coeff_upper[z,m,:])
                # standard_deviation_upper_error[z,m]= standard_deviation_upper [z,m]/pearson_mean_upper[z,m]
              
                # if standard_deviation_upper_error[z,m] <0.01:
                #     pearson_mean_upper[z,m]=pearson_mean_upper[z,m]
                # else:
                #     pearson_mean_upper[z,m]=float('NAN')
                if abs(perfect_linearity_comparison[z,m,i])<0.01:
                     pearson_coeff_upper[z,m,i]=pearson_coeff_upper[z,m,i]
                else:
                     pearson_coeff_upper[z,m,i]=float('NAN')
                     shear_rate_upper[z,m,i]=float('NAN')
                     t_s[z,m,i]=float('NAN')
                     #swap_rate_for_plotting_[z,i,m]='NAN'
                     
# need to get rid of data sets that are no longer complete due to too many NAN points                   
            

#t_s=np.linspace(0,VP_data_upper.shape[3],int(float(no_timesteps)/float(thermo_freq)))   
#%% taking gradient of the shear rate to test for steady state
# for z in range(0,locations_of_non_nan_neg.shape[0]): 
#     for m in range(0,swap_rate.size):
#         for i in range(0,VP_data_upper.shape[3]):
#          shear_rate_grad_upper[z,m,i]=np.gradient(shear_rate_upper[z,m,i],axis=2)
shear_rate_grad_upper=np.gradient(shear_rate_upper,axis=2) # this isnt exact enough need a non-linear curve fit , can then plot the derivative 
z=0
from scipy.optimize import curve_fit
def func(x, a, c):
    return  (a*x)+c               # a * np.exp(-b * x) + c could use an exponential plateau


ydata= shear_rate_upper[z,0,50:]
xdata=t_s[z,0,50:]
popt, pcov = curve_fit(func, xdata,ydata)
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='fit: a=%5.3f, c=%5.3f' % tuple(popt))
plt.plot(xdata,ydata)
plt.show()
#%% plotting to check shear rate is actually steady

f=10
fig=plt.figure(figsize=(13,6))
gs=GridSpec(nrows=1,ncols=2)
fig.suptitle('',size='large',ha='center')

ax1= fig.add_subplot(gs[0,0],projection='3d') #[lowerCorner_x, lowerCorner_y, width, height]
ax2= fig.add_subplot(gs[0,1],projection='3d')

# could I use the mean square error aswell??

plt.rcParams.update(plt.rcParamsDefault)        
   
#t_s=np.linspace(0,VP_data_upper.shape[3],int(float(no_timesteps)/float(thermo_freq)))    

for z in range(0,1):#locations_of_non_nan_neg.shape[0]): 
    for m in range(0,swap_rate.size):
       # for i in range(0,VP_data_upper.shape[3]):
            
            ax1.plot(t_s[z,m,:],swap_rate_for_plotting_[z,:,m],shear_rate_grad_upper[z,m,:],ms=ms)#linestyle='--', marker="x")
            ax1.set_xlabel("$N\Delta t[-]$",fontsize=f,labelpad=f)
            ax1.set_ylabel("$M [-]$", fontsize=f, labelpad=f)
         
            ax1.set_zlabel("$\Delta\dot{\gamma}[-]$", fontsize=f, labelpad=f)
            #ax1.set_title('Upper Velocity Profile')
plt.show()

for z in range(0,1):#locations_of_non_nan_neg.shape[0]): 
    for m in range(0,swap_rate.size):
       # for i in range(0,VP_data_upper.shape[3]):
            
            ax2.plot(t_s[z,m,:],swap_rate_for_plotting_[z,:,m],shear_rate_upper[z,m,:],ms=ms)#linestyle='--', marker="x")
            ax2.set_xlabel("$N\Delta t[-]$",fontsize=f,labelpad=f)
            ax2.set_ylabel("$M [-]$", fontsize=f, labelpad=f)
         
            ax2.set_zlabel("$ \dot{\gamma}[-]$", fontsize=f, labelpad=f)
            #ax1.set_title('Upper Velocity Profile')
plt.show()
#%%
#Lower VP
x= np.array(VP_data_lower_realisation_averaged)
y= np.repeat(np.array([VP_z_data_lower]).T,VP_data_lower.shape[3],1)
pearson_coeff_lower= np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size,VP_data_lower_realisation_averaged.shape[3]))
pearson_mean_lower=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size))    
standard_deviation_lower=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size))                  
perfect_linearity_comparison =np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size))       
standard_deviation_lower_error=np.zeros((locations_of_non_nan_neg.shape[0],swap_rate.size))                           
                              
                              
for z in range(0,locations_of_non_nan_neg.shape[0]): 
    for m in range(0,swap_rate.size):
        for i in range(0,VP_data_lower.shape[3]):
                
                pearson_coeff_lower[z,m,i] =scipy.stats.pearsonr(x[z,m,:,i], y[:,i])[0]
                pearson_mean_lower[z,m]= np.mean(pearson_coeff_lower[z,m,:])
                perfect_linearity_comparison[z,m]=  1-abs(pearson_mean_lower[z,m])
                standard_deviation_lower [z,m]=np.std(pearson_coeff_lower[z,m,:])
                standard_deviation_lower_error[z,m]= standard_deviation_lower [z,m]/pearson_mean_lower[z,m]
                if abs(standard_deviation_lower_error[z,m]) <0.01:
                    pearson_mean_lower[z,m]=pearson_mean_lower[z,m]
                else:
                    pearson_mean_lower[z,m]=float('NAN')
                if perfect_linearity_comparison[z,m]<0.01:
                    pearson_mean_lower[z,m]=pearson_mean_lower[z,m]
                else:
                    pearson_mean_lower[z,m]=float('NAN')
                
#%% plotting pearson vs time to show a steady state is reached
#plt.rcParams['text.usetex'] = True

f=10
fig=plt.figure(figsize=(13,6))
gs=GridSpec(nrows=1,ncols=2)
fig.suptitle('Pearson Correlation coefficent $CC_{p}$ vs $\Delta t_{SRD}/\Delta t_{MD}$ and MPCD bead density $M$',size='large',ha='center')

ax1= fig.add_subplot(gs[0,0],projection='3d') #[lowerCorner_x, lowerCorner_y, width, height]
ax2= fig.add_subplot(gs[0,1],projection='3d')

# could I use the mean square error aswell??

plt.rcParams.update(plt.rcParamsDefault)                 

for z in range(0,locations_of_non_nan_neg.shape[0]):
    for m in range(0,swap_rate.size):
    
        ax1.plot(Solvent_bead_SRD_box_density_cp_1[0,locations_of_non_nan_neg[z,0]],Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]], pearson_mean_upper[z,m],ms=ms,linestyle='--', marker="x")
        ax1.set_xlabel("$M[-]$",fontsize=f,labelpad=f)
        ax1.set_ylabel("$\Delta t_{SRD}/\Delta t_{MD} [-]$", fontsize=f, labelpad=f)
        ax1.set_zlabel("$CC_{p}[-]$", fontsize=f, labelpad=f)
        ax1.set_title('Upper Velocity Profile')
        ax2.plot(Solvent_bead_SRD_box_density_cp_1[0,locations_of_non_nan_neg[z,0]],Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[locations_of_non_nan_neg[z,0],locations_of_non_nan_neg[z,1]], pearson_mean_lower[z,m],ms=ms,linestyle='--', marker="x")        
        ax2.set_xlabel("$M[-]$", fontsize=f, labelpad=f)
        ax2.set_ylabel("$\Delta t_{SRD}/\Delta t_{MD} [-]$", fontsize=f, labelpad=f)
        ax2.set_zlabel("$CC_{p}[-]$", fontsize=f, labelpad=12)
        ax2.set_title('Lower Velocity Profile')
        
        #  ax1.plot( t_s[:],pearson_coeff_upper[z,1,:],Solvent_bead_SRD_box_density_cp_1[0,locations_of_non_nan_neg[z,0]],ms=ms,linestyle='--', marker="+")
       #  ax1.plot( t_s[:],pearson_coeff_upper[z,2,:],Solvent_bead_SRD_box_density_cp_1[0,locations_of_non_nan_neg[z,0]],ms=ms,linestyle='--', marker="o")
       #  ax1.plot( t_s[:],pearson_coeff_upper[z,3,:],Solvent_bead_SRD_box_density_cp_1[0,locations_of_non_nan_neg[z,0]],ms=ms,linestyle='--', marker="v")
       # #  ax1.plot( t_s[:],pearson_coeff_upper[z,2,:],ms=ms,linestyle='--', marker="o")
       #  ax1.plot( t_s[:],pearson_coeff_upper[z,3,:],ms=ms,linestyle='--', marker="v")
       #  ax1.plot( t_s[:],pearson_coeff_upper[z,4,:],ms=ms,linestyle='--', marker="*")
       #  ax1.plot( t_s[:],pearson_coeff_upper[z,5,:],ms=ms,linestyle='--', marker="p")
       #  ax1.plot( t_s[:],pearson_coeff_upper[z,6,:],ms=ms,linestyle='--', marker="h")
       # ax1.set_yscale('log')
        #ax1.grid()
        #ax1.set_yscale('log')
        
#plt.plot(t_s,Mean_pearson_coeff)
plt.tight_layout()
plt.show()    
#%% finding matching inputs np.argwhere(~np.isnan(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg))
locations_of_acceptable_inputs_upper = np.argwhere(~np.isnan(pearson_mean_upper))
locations_of_acceptable_inputs_lower = np.argwhere(~np.isnan(pearson_mean_lower))
valid_tests_M_and_SR =np.zeros((20,2))
#valid_tests_M_and_SR_bool =np.in1d(locations_of_acceptable_inputs_lower,locations_of_acceptable_inputs_upper, assume_unique=True)
#p=np.arange(locations_of_acceptable_inputs_upper.shape[0])[np.in1d(locations_of_acceptable_inputs_upper,locations_of_acceptable_inputs_lower)]        
#valid_tests_M_and_SR=np.reshape(valid_tests_M_and_SR_bool, (locations_of_acceptable_inputs_lower.shape[0],locations_of_acceptable_inputs_lower.shape[1]))
           

for k in range(0,min(locations_of_acceptable_inputs_upper.shape[0],locations_of_acceptable_inputs_lower.shape[0])): 
    for z in range(0,max(locations_of_acceptable_inputs_upper.shape[0],locations_of_acceptable_inputs_lower.shape[0])):
       # for i in range(0,2):
                 
            if max(locations_of_acceptable_inputs_upper.shape[0],locations_of_acceptable_inputs_lower.shape[0]) == locations_of_acceptable_inputs_upper.shape[0]:
                 if  locations_of_acceptable_inputs_lower[k,0]== locations_of_acceptable_inputs_upper[z,0]:
                     
                    if locations_of_acceptable_inputs_lower[k,1]== locations_of_acceptable_inputs_upper[z,1]:
                        
                        valid_tests_M_and_SR[k,0]=np.append(valid_tests_M_and_SR[k,0], locations_of_acceptable_inputs_lower[k,0])[1]
                        valid_tests_M_and_SR[k,1]=np.append(valid_tests_M_and_SR[k,1], locations_of_acceptable_inputs_lower[k,1])[1]
                    else:
                        print()
                    
                     
                 else:
                     print()
                
            else:
                 if locations_of_acceptable_inputs_upper[k,0]== locations_of_acceptable_inputs_lower[z,0]:
                      
                    if locations_of_acceptable_inputs_upper[k,1]== locations_of_acceptable_inputs_lower[z,1]:
                        
                        valid_tests_M_and_SR[k,0]=np.append(valid_tests_M_and_SR[k,0], locations_of_acceptable_inputs_upper[k,0])[1]
                        valid_tests_M_and_SR[k,1]=np.append(valid_tests_M_and_SR[k,1], locations_of_acceptable_inputs_upper[k,1])[1]
                    else:
                        print()
                    
                     
                 else:
                     print()
                    
valid_tests_M_and_SR = valid_tests_M_and_SR[~np.all(valid_tests_M_and_SR== 0, axis=1)]  
# 2nd colum is the swap rate used in this
# can now plug indicies into locations_of_non_nan_neg[] then use those indicies on the original calc
# srd ratio index should be the same for all points 
simulation_point_1_SRD_ratio_indicies =locations_of_non_nan_neg[int(valid_tests_M_and_SR[0,0]),:]
print(simulation_point_1_SRD_ratio_indicies)
simulation_point_1_swap_rate= swap_rate[int(valid_tests_M_and_SR[0,1])]
                                        


#simulation_point_1_SRD_ratio = Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[]
            

#  could plot them with different markers for  each swap rate, and also increase the precision of the axis to spread the data out more 
#%% Code to conserve number of exchnages in a run based on swap freqency
total_number_exchanges=5000 
no_timesteps_after_eq=swap_rate*total_number_exchanges
output_rate_ = swap_rate*100


#%% Simulation Run code with vetted data points 

os.chdir('/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac') 
import mom2numpy
l=5
realisation_index_ =np.linspace(1, l,l)
lamda=mean_free_path_pf_SRD_particles_cp_mthd_1_neg#np.zeros((l,20))
equilibration_timesteps= 1000 # number of steps to do equilibration with 
VP_ave_freq =1000
chunk = 20
no_timesteps=no_timesteps_after_eq
dump_freq=1000 #1000 # if you change the timestep rememebr to chaneg this 
thermo_freq =1000
lamda=mean_free_path_pf_SRD_particles_cp_mthd_1_neg
VP_data_upper=np.zeros((no_timesteps.shape[0],9,no_timesteps.shape[0],5))

for i in range(0,no_timesteps.shape[0]):
               
   VP_data_upper=np.zeros((i,9,int(no_timesteps[i]/thermo_freq),5))
   VP_data_lower=np.zeros((i,9,int(no_timesteps[i]/thermo_freq),5))
  
 #VP_data_upper=np.zeros((9,21,5))





#VP_data_lower=np.zeros((no_timesteps.shape[0],9,int(no_timesteps/thermo_freq),5))
#%%
log_files=np.zeros((valid_tests_M_and_SR.shape[0],number_of_realisations,int((no_timesteps+equilibration_timesteps)/thermo_freq),4))
#VP_data_lower=np.zeros((9,21,5))
#log_file = np.array([[[]]])

simulation_file="Simulation_run_folder"
number_of_realisations=3
log_files=np.zeros((valid_tests_M_and_SR.shape[0],number_of_realisations,int((no_timesteps+equilibration_timesteps)/thermo_freq),4))
momentum_data =np.zeros((valid_tests_M_and_SR.shape[0],number_of_realisations))
VP_data_upper=np.zeros((valid_tests_M_and_SR.shape[0],number_of_realisations,9,int(no_timesteps/thermo_freq)))
VP_data_lower=np.zeros((valid_tests_M_and_SR.shape[0],number_of_realisations,9,int(no_timesteps/thermo_freq)))
z=0
print(locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),:])
l=8 # swap number index

#%% Simulation Run loop 
for z in range(0,1):#valid_tests_M_and_SR.shape[0]):
    for j in range(0,1):#number_of_realisations):#5):# for now just use one realisation 
       # for m in range(0,1):#swap_rate.size):#range(0,1):       #range(0,swap_rate.size):
          #  for k in range(0,1):#,swap_number.size):
                        
                        
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]]))
                        print("Number of SRD particles:",no_SRD)
                        
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]])
                        print("Scaled mass of an SRD:",mass_SRD)
                        box_size = str(box_side_length_scaled)
                        print("Scaled Box side length:",box_size) 
                        timestep_input= str(scaled_timestep)
                        print("Scaled timestep:",timestep_input)
                        chunk = 20 # number of chunks to use for VP averaging
                        #if str(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i])=='nan' :
                           # continue 
                        #else:
                            
                        SRD_MD_ratio=str(int(np.round( Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]])))
                        
                        #    SRD_MD_ratio=str(int(np.ceil( Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i])))
                             
                        # SRD_MD_ratio=str(int(np.ceil( Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i])))
                        print("Number of MD steps per SRD step:",SRD_MD_ratio) # select one input 
                        lamda[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]] = mean_free_path_pf_SRD_particles_cp_mthd_1_neg[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]]
                        grid_size = str(SRD_box_size_wrt_solid_beads[0,locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]])
                        print("SRD grid size:",grid_size)
                        #print("alpha:"+str(alpha[0,i]))
                        shear_rate =str(scaled_gamma_dot[0,1]) # which one do i use ?? or the average 
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps_ = str(no_timesteps[int(valid_tests_M_and_SR[z,1])])
                        #no_timesteps_ =float(no_timesteps)
                        temp_=str(scaled_temp)# arbit
                        
                        #t_damp = str(t_damp_SRD[j])
                        r_part_sim = str(r_particle_scaled)
                        print("Scaled temp:",temp_)
                         # needs to be calculated to run a specific amount of strain 
                       # total_strain = float(shear_rate) * scaled_timestep * no_timesteps_ 
                       # print("Total strain:",total_strain)
                       # print("Expected dynamic viscosity (dimensionless):",scaled_viscosity)
                        rand_int =str(np.random.randint(0, 10000))
                        rand_int_1 =str( np.random.randint(0, 10000))
                        rand_int_2 =str(np.random.randint(0, 10000))
                        
                   
                
                # Run code 
                
                        os.chdir('/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder') 
                        # /Volumes/Backup\ Plus/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/Simulation\ run\ folder\  
                        #/Volumes/Backup\ Plus/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/Simulation\ run\ folder\
                       
                    
                    
                   
                        os.system('mpirun -np 8 ./lmp_mpi -var output_rate_ '+str(output_rate_[int(valid_tests_M_and_SR[z,1])])+' -var swap_rate '+str(swap_rate[int(valid_tests_M_and_SR[z,1])])+' -var swap_number '+str(swap_number[l])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var r_part_sim '+r_part_sim+' -var temp_ '+temp_+' -var lambda '+str(lamda[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]])+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var shear_rate ' +shear_rate+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps_+' -in no_wall_pure_SRD_sim_var_inputs_td_var_no_tstat_no_rescale.file')
    
        
                 #   os.system('mpirun -np 8 ./lmp_mpi -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var r_part_sim '+r_part_sim+' -var temp_ '+temp_+' -var lambda '+str(lamda[j,i])+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var shear_rate ' +shear_rate+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -in no_wall_pure_SRD_sim_var_inputs_td_var_fix_deform_version.file')
                        Path_2_VP="/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/"+simulation_file
                       
                        #thermo_vars ="         KinEng          Temp          TotEng    "
                        
                        VP_output_col_count = 4 
                        

                        realisation_name="vel.profile_no_tstat__no_rescale_"+str(rand_int)+"_"+str(realisation_index_[j])+"_"+str(no_SRD)+"_"+str(box_size)+"_"+str(timestep_input)+"_"+str(SRD_MD_ratio)+"_"+str(dump_freq)+"_"+str(thermo_freq)+"_"+str(no_timesteps)+"_T_"+str(temp_)+"_lbda_"+str(lamda[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]])+"_SR_"+str(swap_rate[int(valid_tests_M_and_SR[z,1])])+"_SN_"+str(swap_number[l])  
                                                                           #vel.profile_no_tstat__no_rescale_9996_1.0_7507_16.0_0.001_381_1000_1000_20000_T_5e-05_lbda_0.00934921866630143_SR_3_SN_1
                        #realisation_name_="log.no_wall_no_tstat_no_rescale_"+str(rand_int)+"_"+str(realisation_index_[j])+"_"+str(no_SRD)+"_"+str(box_size)+"_"+str(timestep_input)+"_"+str(SRD_MD_ratio)+"_"+str(dump_freq)+"_"+str(thermo_freq)+"_"+str(no_timesteps)+"_T_"+str(temp_)+"_lbda_"+str(lamda[simulation_point_1_SRD_ratio_indicies[0],simulation_point_1_SRD_ratio_indicies[1]])+"_SR_"+str(simulation_point_1_swap_rate)+"_SN_"+str(swap_number[0])  
                        #Path_2_dump="/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/"+simulation_file
                       
                        VP_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[0]
                        VP_z_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[1]
                       
                        VP_data_upper[z,j,:,:] = VP_data[1:10,:]
                        VP_z_data_upper = VP_z_data[1:10].astype('float64')* box_side_length_scaled
                        VP_data_lower[z,j,:,:] = VP_data[11:,:]
                 
                        VP_z_data_lower = VP_z_data[11:].astype('float64') * box_side_length_scaled
                        
                        Path_2_log=Path_2_VP
                        realisation_name="log.no_wall_no_tstat_no_rescale_"+str(rand_int)+"_"+str(realisation_index_[j])+"_"+str(no_SRD)+"_"+str(box_size)+"_"+str(timestep_input)+"_"+str(SRD_MD_ratio)+"_"+shear_rate+"_"+str(dump_freq)+"_"+str(thermo_freq)+"_"+str(no_timesteps)+"_T_"+str(temp_)+"_lbda_"+str(lamda[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]])+"_SR_"+str(swap_rate[int(valid_tests_M_and_SR[z,1])])+"_SN_"+str(swap_number[l])  
                        #log.no_wall_no_tstat_no_rescale_742_1.0_7507_16.0_0.001_381_0.01_1000_1000_20000_T_5e-05_lbda_0.00934921866630143_SR_3_SN_1
                        thermo_vars ="         KinEng          Temp          TotEng    "
                        log_files[z,j,:,:] = log2numpy.log2numpy(Path_2_log,thermo_vars,realisation_name)[0]
                        # fix the realisation index line check it works properly, might not be calling the correct name 
                        realisation_name="mom._no_tstat__no_rescale_"+str(rand_int)+"_"+str(realisation_index_[j])+"_"+str(no_SRD)+"_"+str(box_size)+"_"+str(timestep_input)+"_"+str(SRD_MD_ratio)+"_"+str(dump_freq)+"_"+str(thermo_freq)+"_"+str(no_timesteps)+"_T_"+str(temp_)+"_lbda_"+str(lamda[locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),0],locations_of_non_nan_neg[int(valid_tests_M_and_SR[z,0]),1]])+"_SR_"+str(swap_rate[int(valid_tests_M_and_SR[z,1])])+"_SN_"+str(swap_number[l])
                        #mom._no_tstat__no_rescale_9671_1.0_7507_16.0_0.001_381_1000_1000_20000_T_5e-05_lbda_0.00934921866630143_SR_3_SN_1
                        momentum_data[z,j] = mom2numpy.mom2numpy(realisation_name,Path_2_log)
#%%Processing data 

"""
This section will now average the VP data, check the pearson coeff is approx 1, then extract a shear rate, the  calculate the flux
then finally get the viscosity 
"""

VP_data_lower_realisation_averaged = np.mean(VP_data_lower,axis=1) 
# 4 indicates to take the mean of the 
VP_data_upper_realisation_averaged = np.mean(VP_data_upper,axis=1) 
x= np.array(VP_data_lower_realisation_averaged)
x_=np.array(VP_data_upper_realisation_averaged)
y= np.repeat(np.array([VP_z_data_lower]).T,VP_data_lower.shape[3],1)
y_=np.repeat(np.array([VP_z_data_upper]).T,VP_data_upper.shape[3],1)
pearson_coeff_lower_post_run= np.zeros((valid_tests_M_and_SR.shape[0],VP_data_lower_realisation_averaged.shape[2]))
lower_data_slope=np.zeros((valid_tests_M_and_SR.shape[0],VP_data_lower_realisation_averaged.shape[2]))
pearson_coeff_upper_post_run= np.zeros((valid_tests_M_and_SR.shape[0],VP_data_upper_realisation_averaged.shape[2]))
upper_data_slope= np.zeros((valid_tests_M_and_SR.shape[0],VP_data_upper_realisation_averaged.shape[2]))

for z in range(0,valid_tests_M_and_SR.shape[0]):
       #for j in range(0,number_of_realisations):
          for i in range(0,VP_data_lower.shape[3]):
              
                 pearson_coeff_lower_post_run[z,i]=scipy.stats.pearsonr(x[z,:,i], y[:,i])[0]
                 lower_data_slope[z,i]=scipy.stats.linregress(y[:,i],x[z,:,i] ).slope
                 pearson_coeff_upper_post_run[z,i]=scipy.stats.pearsonr(x_[z,:,i], y_[:,i])[0]
                 upper_data_slope[z,i]=scipy.stats.linregress( y_[:,i],x_[z,:,i]).slope



t_s=np.linspace(0,VP_data_upper.shape[3],int(float(no_timesteps)/float(thermo_freq))).T    

#%%
plt.title("Test Run showing the shear rate reaching a steady state")
plt.xlabel("$N\Delta t [-]$")
plt.ylabel("$\dot{\gamma}[-]$")
plt.plot(t_s,upper_data_slope[0,:])

plt.show()

























#%% old velocity profile code
 
    #                 Path_2_VP="/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/"+simulation_file
    #                 VP_output_col_count = 4 
    #                 realisation_name="vel.profile_"+str(rand_int)+"_"+str(realisation_index_[j])+"_"+str(no_SRD)+"_"+str(box_size)+"_"+str(timestep_input)+"_"+str(SRD_MD_ratio)+"_"+str(dump_freq)+"_"+str(thermo_freq)+"_"+str(no_timesteps)+"_T_"+str(temp_)+"_lbda_"+str(lamda[z,i])+"_SR_"+str(swap_rate[m])+"_SN_"+str(swap_number[k])  
    #                 #realisation_name="vel.profile_3302_1.0_72000_1.0_0.001_5_1000_1000_200000_T_1.0_lbda_3.900861786968397e-05_SR_3_SN_50"
    #                 VP_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[0]
    #                 VP_z_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[1]
    # #          # split the profile into mirror images  
    #                 VP_data_upper = VP_data[1:10,:]
    #                 VP_z_data_upper = VP_z_data[1:10].astype('float64')* box_side_length_scaled
    #                 VP_data_lower = VP_data[11:,:]
    #                 VP_z_data_lower = VP_z_data[11:].astype('float64') * box_side_length_scaled
    #                 # upper profile 
    #          # looking at the time series of velocity profiles 
    #                 # for i in range(0,VP_data_upper.shape[1]):
    #                 #     x=VP_z_data_upper 
    #                 #     y=VP_data_upper[:,i]
    #                 #     #plt.suptitle('Upper Velocity Profile(SR:'+str(swap_rate)+',SN:'+str(swap_number)+') $\delta t_{MD}='+str(scaled_timestep)+'$',size='large', wrap=True )
    #                 #    # plt.xlabel('$\frac{z}{l^{*}} [-]$', size='large')
    #                 #     #plt.ylabel( '$\nu_{x}^{*} [-]$', rotation='horizontal',ha='right',size='large')
    #                 #     plt.plot(x,y)
    #                 #     plt.show()
    #                 # # lower profile 
    #                 # for i in range(0,VP_data_upper.shape[1]):
    #                 #     x=VP_z_data_lower 
    #                 #     y=VP_data_lower[:,i]
    #                 #     #plt.suptitle('Lower Velocity Profile(SR:'+str(swap_rate)+',SN:'+str(swap_number)+') $\delta t_{MD}='+str(scaled_timestep)+'$',size='large', wrap=True )
    #                 #     #plt.set_xlabel('$\frac{z}{l^{*}} [-]$', size='large')
    #                 #     #plt.set_ylabel( '$\nu_{x}^{*} [-]$', rotation='horizontal',ha='right',size='large')
    #                 #     plt.plot(x,y)
    #                 #     plt.show()
    #                 # need to take the velocity profile where the total energy of the system has reached a steady state, can then take an average value for each chunk pover the steady state
    #                 # first i need to fix the log file reader so it works with LAMMPS Jun 2022
    #                 # using velocity profile data to check the shear rate 
    #                 # only using last 5 runs for now and cutting off the first piece of data in each as there is discontinuity over the middle 
    #                 number_of_VP_outputs = int(float(no_timesteps)/float(VP_ave_freq))
    #                 upper_y_1= VP_data_upper[8,number_of_VP_outputs-5:]
    #                 upper_y_2= VP_data_upper[0,number_of_VP_outputs-5:]
    #                 upper_z_data_x_1= float(VP_z_data_upper[8])
    #                 upper_z_data_x_2= float(VP_z_data_upper[0])
    #                 delta_z_upper =upper_z_data_x_2-upper_z_data_x_1
    #                 upper_calculated_shear_rate = (upper_y_2-upper_y_1)/(upper_z_data_x_2-upper_z_data_x_1)
                    
    #                 lower_y_1= VP_data_lower[8,number_of_VP_outputs-5:]
    #                 lower_y_2= VP_data_lower[0,number_of_VP_outputs-5:]
    #                 lower_z_data_x_1= float(VP_z_data_lower[8])
    #                 lower_z_data_x_2= float(VP_z_data_lower[0])
    #                 delta_z_lower =lower_z_data_x_2-lower_z_data_x_1
    #                 lower_calculated_shear_rate = (lower_y_2-lower_y_1)/(lower_z_data_x_2-lower_z_data_x_1)
                    
    #                 #m,k
 
    #                 shear_rate_average[m,i,] = (np.mean(abs(upper_calculated_shear_rate))+np.mean(abs(lower_calculated_shear_rate)) )*0.5        
    #                 #plotting upper and lower profile 
                   
#%% loading individual VPs
#chunk = 20 # number of chunks to use for VP averaging 
#equilibration_timesteps= 1000 # number of steps to do equilibration with 
#VP_ave_freq =1000
#no_SRD=str(72000)
# import velP2numpy
# #no_timesteps = str(200000)
# simulation_file="Simulation_run_folder"
# chunk = 20 # number of chunks to use for VP averaging 
# equilibration_timesteps= 1000 # number of steps to do equilibration with 
# VP_ave_freq =1000
# no_SRD=str(28000)
# no_timesteps = str(200000)
# Path_2_VP="/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/"+simulation_file
# VP_output_col_count = 4 
# realisation_name="vel.profile_8945_1.0_28000_16.0_0.001_48_1000_1000_200000_T_5e-05_lbda_0.0022412863997712524_SR_3_SN_50"

# realisation_name="vel.profile_8945_1.0_28000_16.0_0.001_48_1000_1000_200000_T_5e-05_lbda_0.0022412863997712524_SR_3_SN_50"
# VP_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[0]
# VP_z_data = velP2numpy.velP2numpy(Path_2_VP,chunk,realisation_name,equilibration_timesteps,VP_ave_freq,no_SRD,no_timesteps,VP_output_col_count)[1]
# #          # split the profile into mirror images  
#%% Calculating shear profiles
# VP_data_upper = VP_data[:10,:]
# VP_z_data_upper = VP_z_data[:10].astype('float64')* box_side_length_scaled
# VP_data_lower = VP_data[10:,:]
# VP_z_data_lower = VP_z_data[10:].astype('float64') * box_side_length_scaled
 
# number_of_VP_outputs = int(float(no_timesteps)/float(VP_ave_freq))
# upper_y_1= VP_data_upper[8,number_of_VP_outputs-5:]
# upper_y_2= VP_data_upper[0,number_of_VP_outputs-5:]
# upper_z_data_x_1= float(VP_z_data_upper[8])
# upper_z_data_x_2= float(VP_z_data_upper[0])
# delta_z_upper =upper_z_data_x_2-upper_z_data_x_1
# upper_calculated_shear_rate = (upper_y_2-upper_y_1)/(upper_z_data_x_2-upper_z_data_x_1)

# lower_y_1= VP_data_lower[8,number_of_VP_outputs-5:]
# lower_y_2= VP_data_lower[0,number_of_VP_outputs-5:]
# lower_z_data_x_1= float(VP_z_data_lower[8])
# lower_z_data_x_2= float(VP_z_data_lower[0])
# delta_z_lower =lower_z_data_x_2-lower_z_data_x_1
# lower_calculated_shear_rate = (lower_y_2-lower_y_1)/(lower_z_data_x_2-lower_z_data_x_1)

# #m,k
# shear_rate_average = (np.mean(abs(upper_calculated_shear_rate))+np.mean(abs(lower_calculated_shear_rate)) )*0.5      

#%% Velocity profiles test
#upper
for i in range(0,VP_data_upper.shape[1]):
    x=VP_z_data_upper 
    y=VP_data_upper[:,i]
    #plt.suptitle('Upper Velocity Profile(SR:'+str(swap_rate)+',SN:'+str(swap_number)+') $\delta t_{MD}='+str(scaled_timestep)+'$',size='large', wrap=True )
    # plt.xlabel('$\frac{z}{l^{*}} [-]$', size='large')
    #plt.ylabel( '$\nu_{x}^{*} [-]$', rotation='horizontal',ha='right',size='large')
    plt.plot(x,y)
    plt.show()
    
#%% upper new code
x_data_cp_1 = np.linspace(1000,201000,201)
#print(x_data_cp_1)
x_temporary=np.array([[VP_z_data_upper]])

y_data_cp_1 = np.repeat(x_temporary,x_data_cp_1.size,axis=1)
#print(y_data_cp_1)
z_data_cp_1_pos = VP_data_upper.T
#z_data_cp_1_neg = mean_free_path_pf_SRD_particles_cp_mthd_1_neg
#%% new vp code for upper profile 
fig=plt.figure(figsize=(13,6))
gs=GridSpec(nrows=1,ncols=1)
# ax1= fig.add_subplot(1,1,1,adjustable='box', aspect=1)
# ax2= fig.add_subplot(1,2,2,adjustable='box', aspect=1)
#fig, (ax1,ax2)=plt.subplots(1,2)
file="Upper Velocity Profile for "+realisation_name
fig.suptitle('Upper Velocity Profile for '+realisation_name ,size='large',ha='center')

ax1= fig.add_subplot(gs[0]) #[lowerCorner_x, lowerCorner_y, width, height]
#ax2= fig.add_subplot(gs[1])     #[lowerCorner_x, lowerCorner_y, width, height]

# if you need to reset stle as plt.style is not responding use the line below
plt.rcParams.update(plt.rcParamsDefault) 
#seaborn-notebook
#plt.style.use('seaborn-paper')
#plt.subplots_adjust(hspace=5)
for i in range (0,x_data_cp_1.size-1):
    #for j in range (0,4): 
     ax1.plot( y_data_cp_1[0,i,:], z_data_cp_1_pos[i,:],ms=ms,linestyle='--', marker="x")#label='M+={}'.format(x_data_cp_1[i]),
     #ax1.plot( y_data_cp_1[i,:], z_data_cp_1_neg[i,:],ms=ms,label='M-={}'.format(x_data_cp_1[i]),linestyle='--', marker="+")
    # ax1.legend(x_data_cp_1)
     ax1.set_xscale('linear')
     ax1.set_yscale('linear')
     #ax1.set_title('SRD Coupling Method 1', size='x-large')
     ax1.set_xlabel('$ L^{*} [-]$',size='large')
     ax1.set_ylabel( '$u^{*} [-]$', rotation='horizontal',ha='right',size='large')
    # ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left',borderaxespad=1)
     
     ax1.grid()
    # ax1.imshow( aspect='auto')
#on y asix
     

# =============================================================================
# for i in range (0,number_of_M_cp_2):
#     #for j in range (0,4): 
#      ax2.plot( y_data_cp_2[i,:], z_data_cp_2[i,:],label='M={}'.format(x_data_cp_2[i]),linestyle='--', marker='o')
#      ax2.legend(x_data_cp_1)
#    
#      ax2.set_xscale('log')
#      ax2.grid()
#      ax2.set_yscale('log')
#      ax2.set_title('SRD Coupling Method 2', size='x-large') 
#      ax2.set_xlabel('$\lambda/a$',size='x-large')
#      ax2.set_ylabel( '$\delta t_{SRD}/\delta t_{MD}$', rotation='horizontal',ha='right',size='x-large')
#      ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=1)
#      #ax2.set_aspect(aspect='equal',adjustable='box')
# 
#      
# =============================================================================
plt.tight_layout()

path= os.chdir("/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder/dumps to show helen ")
#plt.savefig(file+'.png')
plt.show()
#%% upper new code
x_data_cp_1 = np.linspace(1000,201000,201)
#print(x_data_cp_1)
x_temporary=np.array([[VP_z_data_lower]])

y_data_cp_1 = np.repeat(x_temporary,x_data_cp_1.size,axis=1)
#print(y_data_cp_1)
z_data_cp_1_pos = VP_data_lower.T
#%% new vp code for lower profile 
fig=plt.figure(figsize=(13,6))
gs=GridSpec(nrows=1,ncols=1)
# ax1= fig.add_subplot(1,1,1,adjustable='box', aspect=1)
# ax2= fig.add_subplot(1,2,2,adjustable='box', aspect=1)
#fig, (ax1,ax2)=plt.subplots(1,2)
file="Lower Velocity Profile for "+realisation_name
fig.suptitle('Lower Velocity Profile for '+realisation_name ,size='large',ha='center')

ax1= fig.add_subplot(gs[0]) #[lowerCorner_x, lowerCorner_y, width, height]
#ax2= fig.add_subplot(gs[1])     #[lowerCorner_x, lowerCorner_y, width, height]

# if you need to reset stle as plt.style is not responding use the line below
plt.rcParams.update(plt.rcParamsDefault) 
#seaborn-notebook
#plt.style.use('seaborn-paper')
#plt.subplots_adjust(hspace=5)
for i in range (0,x_data_cp_1.size-1):
    #for j in range (0,4): 
     ax1.plot( y_data_cp_1[0,i,:], z_data_cp_1_pos[i,:],ms=ms,linestyle='--', marker="x")#label='M+={}'.format(x_data_cp_1[i]),
     #ax1.plot( y_data_cp_1[i,:], z_data_cp_1_neg[i,:],ms=ms,label='M-={}'.format(x_data_cp_1[i]),linestyle='--', marker="+")
    # ax1.legend(x_data_cp_1)
     ax1.set_xscale('linear')
     ax1.set_yscale('linear')
     #ax1.set_title('SRD Coupling Method 1', size='x-large')
     ax1.set_xlabel('$ L^{*} [-]$',size='large')
     ax1.set_ylabel( '$u^{*} [-]$', rotation='horizontal',ha='right',size='large')
    # ax1.legend(bbox_to_anchor=(1.02, 1), loc='Lower left',borderaxespad=1)
     
     ax1.grid()
    # ax1.imshow( aspect='auto')
#on y asix
     

# =============================================================================
# for i in range (0,number_of_M_cp_2):
#     #for j in range (0,4): 
#      ax2.plot( y_data_cp_2[i,:], z_data_cp_2[i,:],label='M={}'.format(x_data_cp_2[i]),linestyle='--', marker='o')
#      ax2.legend(x_data_cp_1)
#    
#      ax2.set_xscale('log')
#      ax2.grid()
#      ax2.set_yscale('log')
#      ax2.set_title('SRD Coupling Method 2', size='x-large') 
#      ax2.set_xlabel('$\lambda/a$',size='x-large')
#      ax2.set_ylabel( '$\delta t_{SRD}/\delta t_{MD}$', rotation='horizontal',ha='right',size='x-large')
#      ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=1)
#      #ax2.set_aspect(aspect='equal',adjustable='box')
# 
#      
# =============================================================================
plt.tight_layout()

path= os.chdir("/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder/dumps to show helen ")
#plt.savefig(file+'.png')
plt.show()
#%% Velocity profiles test
# # lower profile 
for i in range(0,VP_data_upper.shape[1]):
    x=VP_z_data_lower 
    y=VP_data_lower[:,i]
    #plt.suptitle('Lower Velocity Profile(SR:'+str(swap_rate)+',SN:'+str(swap_number)+') $\delta t_{MD}='+str(scaled_timestep)+'$',size='large', wrap=True )
    #plt.set_xlabel('$\frac{z}{l^{*}} [-]$', size='large')
    #plt.set_ylabel( '$\nu_{x}^{*} [-]$', rotation='horizontal',ha='right',size='large')
    plt.plot(x,y)
    plt.show()


#%% Checking energy is steady state    

#Checking energy is at steady state my by reading in all the data files 
import log2numpy


            #r_part_sim=str(5.0)
simulation_file="Simulation_run_folder"
# realisation line below could be made to contain vars 
realisation_name = "log.no_wall_9153_1.0_28000_16.0_0.001_48_0.01_1000_1000_200000_T_5e-05_lbda_0.0022412863997712524_SR_3_SN_1"
#realisation_name = "log.no_wall_rand_int_"+str(rand_int)+"_RI_"+str(realisation_index_[j])+"_lambda_"+str(lamda[j,i])+"_"+no_SRD+"_"+box_size+"_"+timestep_input+"_"+SRD_MD_ratio+"_"+shear_rate+"_"+dump_freq+"_"+thermo_freq+"_"+no_timesteps+"_T_"+temp_+"_"+r_part_sim
                                      #log.fixed_wall_t_27706_16.0_1e-05_34_1.0203153383802774_1000_100_300000_T_1.0_td_0.1_1
thermo_vars ="         KinEng          Temp          TotEng"
 
log_file=[]
Path_2_dump="/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/"+simulation_file
Path_2_log=Path_2_dump


log_test = log2numpy.log2numpy(Path_2_log,thermo_vars,realisation_name)[0]
SRD_temp_lambda =log2numpy.log2numpy(Path_2_log,thermo_vars,realisation_name)[1]
SRD_temp_lambda = SRD_temp_lambda.split()
#SRD_temp[j,i] = SRD_temp_lambda[5] # 5 should always be the temp index 

# stress vs timestep cp_1

#%%
# =============================================================================
x = log_test[:,0]
y=log_test[:,3]
plt.figure(figsize=(13,6))
plt.suptitle('Total Energy vs Timestep '+realisation_name,wrap=True) #(alpha='+str(alpha[0,i])+','+no_SRD+' SRDs,'+str(round(total_strain))+' strain units, timestep='+timestep_input+ '$t^{-1}$),temp='+temp_+',tdamp='+t_damp+', RI:'+str(realisation_index_[j]),wrap=True )
plt.xlabel(' $\Delta t_{MD}[-]$')
plt.ylabel('$E_{total}[-]$',rotation='horizontal')

plt.tight_layout
plt.plot(x,y)
plt.show()

#plt.savefig('Total Energy vs Timestep '+realisation_name+'.png')#'fw_rand_int_'+str(rand_int)+'_RI_'+str(realisation_index_[j])+' alpha_'+str(alpha[0,i])+'_'+no_SRD+'_SRDs_SRD_MD_ratio_'+SRD_MD_ratio+'_'+str(round(total_strain))+'_'+shear_rate+'_strain_tstep_'+timestep_input+ 'cp_1_T_'+temp_+'r_particle_'+r_part_sim+'.png')
# plt.close('all')
    # determine average stress after steady state 

   
#realisation_name = "log.fixed_wall_"+str(realisation_index_[j])+"_"+no_SRD+"_"+box_size+"_"+timestep_input+"_"+SRD_MD_ratio+"_"+shear_rate+"_"+dump_freq+"_"+thermo_freq+"_"+no_timesteps+"_T_"+temp_+"_"+r_part_sim
# log_data = log2numpy.log2numpy(Path_2_log,thermo_vars,realisation_name)[0]
# SS_cutoff =0 #100000 from plot # input("Please input Steady state cut off")
# stress_raw=log_data[SS_cutoff:,1]
# stress_ave[j,i] = np.mean(stress_raw)
# print(stress_ave) 
# # not 100% sure on the method below
# apparent_visc[j,i] = (abs(stress_ave[j,i]) /float(shear_rate)) * SRD_mass_scale_parameter/(lengthscale_parameter* timescale_parameter ) 
# print(apparent_visc[j,i])

#%% 
# x = t_damp_SRD[:]
# y= stress_ave[:,i]
# plt.suptitle('Apparent Stress vs Temperature damping frequency(alpha='+str(alpha[0,i])+','+no_SRD+' SRDs, timestep='+timestep_input+ '$t^{-1}$',wrap=True )
# plt.ylabel('$\sigma$ ''$[ML^{-2}t^{-1}]$')
# plt.xlabel('$ tdamp [-] $')

# plt.plot(x,y)
# plt.show()
# #%%
# nd_visc = 80.359873/100000# 0.66295093#example visc from previous outputs 
# dimensional_visc = energy_parameter * timescale_parameter * nd_visc / lengthscale_parameter**(3)
# print(dimensional_visc) # need to check this was in pascal s
# nd_lambda = 3.575e-05 * lengthscale_parameter
# print(nd_lambda)
# nd_lambda = 1.7875e-09
# lamda_d  = SRD_box_size_wrt_solid_beads * nd_lambda
# print(lamda)
# #%%
# apparent_visc_test = (abs(-0.8) /float(1)) * SRD_mass_scale_parameter/(lengthscale_parameter* timescale_parameter ) 
# print(apparent_visc_test)
# # -0.8 was -7.5e-5
# #%%

# alpha=np.array([np.linspace(0.00006435,0.00007865,15)])
# lamda_ = scaled_SRD_grid_size * alpha

# T_srd = ((lamda_)**2) * SRD_mass_with_solid_cp_mthd_1_scaled[0,1] / (scaled_SRD_timestep_with_solid_cp_mthd_1[1,:])**(2) *scaled_k_b
# #%%
# x = lamda[0,:]
# y=SRD_temp[0,:]
# y_=T_srd[0,:]
# plt.suptitle('$T^{*}_{SRD}[-]$ vs $\lambda^{*}[-]$ ,'+no_SRD+' SRDs, timestep='+timestep_input+ '$t^{-1}$), RI:'+str(realisation_index_[j]),wrap=True )
# plt.xlabel('$\lambda^{*}[-]$')
# plt.ylabel('$T^{*}_{SRD}[-]$')
# plt.legend()
# plt.plot(x,y)
# plt.plot(x,y_)
# plt.show()

#%% using the dump file reader to output the angular velocity vector 

# import dump2numpy

# dump_start_line='ITEM: ATOMS id type x y z vx vy vz'
# filename=''
# dump_realisation_name= 'no_wall_1.0_114688_16.0_0.0001_120_1.0_1_1000_100000_T_1.0_1.0_td_5.0.dump'
# dump_data = dump2numpy.dump2numpy(dump_start_line, Path_2_dump, simulation_file, filename, dump_realisation_name)

# in this file we are looking for particle 1, 19451, 38539, 57655, 76853,95815
#%% omega vector 
# this code is broken , need to find the index opf the row where the solid particles data is rethink this 
# omega_vector = np.array([])
# solid_pattern = re.compile(b"1")
# for i in range(0,np.shape(dump_data)[0]):
#     omega_vector_search = dump_data[i][:,1] 
#     omega_vector_row_index = np.where(omega_vector_search==1)



