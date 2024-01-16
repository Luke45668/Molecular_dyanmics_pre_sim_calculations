#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday 25th Sept 2023

This script will does all the pre calcs for the SRD fluid mappings and then produces the simulation run files for myriad. 


@author: lukedebono
"""
#%%
#from msilib import MSIMODIFY_INSERT_TEMPORARY
import os
import numpy as np

import matplotlib.pyplot as plt
import regex as re
import pandas as pd

# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Helvetica"
# })
from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
#import seaborn as sns
import math as m
import scipy.stats
from datetime import datetime
from sim_file_producer_SRD import *

rho_solid=1200
#%% H20 parameters
fluid_name='H20'
real_fluid_density = 1000 # H20 
T_K=300 # fluid temp at density listed

#%% Argon 
fluid_name='Ar'
real_fluid_density = 1426.9 # H20 
T_K=86.5 # fluid temp at density listed
#%% Nitrogen 
fluid_name='N2'
real_fluid_density =  847 # H20 
T_K=72.2  # fluid temp at density listed
#%% hexane 
fluid_name='C6H14'
real_fluid_density = 700# H20 
T_K=311 # fluid temp at density listed
#%% cyclohexane 
fluid_name='C6H12'
real_fluid_density = 764.95 # H20 
T_K=307.65 # fluid temp at density listed
#%% calcuating lengthscale coefficient 
n_particle=2
phi=np.array([0.05,0.03,0.01,0.008,0.005,0.003,0.002,0.0015,0.0013,0.0012])
fina_size=300 # need to increase this if you want bigger boxes
box_size_bar=np.linspace(1,fina_size,fina_size)

def func_phi(a,n_particle,L,phi):

  return  (4/3)*(np.pi*n_particle/(phi*((L*a)**3))) - 1

#NOTE: for adding the solid particle this should actually be 0.5 x radius , we did 0.25 x radius 
a=np.array([0.5]) # making this bigger reduces the total number of boxes and the particle count

res=np.zeros((phi.size,fina_size,a.size))

for k in range(0,phi.size):
    for i in range(0,a.size):
        for j in range(0,fina_size):
            
            res[k,j,i]=func_phi(a[i],n_particle,box_size_bar[j],phi[k])
         
sign_change_array=np.zeros((phi.size,2,a.size))

for k in range(0,phi.size):
    for i in range(0,a.size):
        for j in range(0,fina_size-1):

    # checking for successive opposite index
            if  res[k,j,i] > 0 and res[k,j+1,i] < 0 or res[k,j,i] < 0 and res[k,j+1,i] > 0:
                sign_change_array[k,0,i]= j
                sign_change_array[k,1,i]= j+1
            else:
                 continue


final_sign_change_array=np.zeros((phi.size,1,a.size))
for k in range(0,phi.size):
    for i in range(0,a.size):
     if  np.abs(res[k,int(sign_change_array[k,0,i]),i]) < np.abs(res[k,int(sign_change_array[k,1,i]),i]):
        
        final_sign_change_array[k,0,i]=int(sign_change_array[k,0,i])
     else:
        final_sign_change_array[k,0,i]=int(sign_change_array[k,1,i])
    
integer_box_size= final_sign_change_array.flatten()
         
k_for_each_phi= np.cbrt((4/3)*np.pi*n_particle/(phi * integer_box_size**3))
 
# check 
def func_phi_check(n_particle,k_for_each_phi,integer_box_size):
    
    return (4/3)* np.pi * n_particle / ((k_for_each_phi*integer_box_size)**3)

phi_check=func_phi_check(n_particle,k_for_each_phi,integer_box_size)

if phi.all()==phi_check.all():
    print("Box sizing consistent")
else: 
    print("check calculation")

# input calc 
k_b=1.38e-23
collision_cell_size_bar= 1
fluid_particle_mass_bar= 1 
simulation_temp_bar=1
fluid_particle_number_density= 5 
energy_scale= T_K * k_b / simulation_temp_bar
r_particle=2.5e-5
lengthscale= k_for_each_phi * r_particle
mass_scale= real_fluid_density * (lengthscale**3) /  fluid_particle_number_density
time_scale= np.sqrt((mass_scale*(lengthscale**2))/(T_K*k_b))
#nu_bar=0.52
#nu_bar=0.7
nu_bar=0.9
#
# need to change this if you change fluid particle number density 
M=fluid_particle_number_density
def gamma_1(M):
    gamma_1=1-((1-np.exp(-M))/M)
    return gamma_1
def gamma_2(M):
    gamma_2= (M+2)/(M-1)
    return gamma_2



gamma_1= gamma_1(M)
gamma_2= gamma_2(M)
gamma_1_2=gamma_1* gamma_2 *(1/18)

#plotting timestep vs nu_bar
def collision_time(nu_bar,gamma_1,gamma_2):
    collision_time_negative_bar= (nu_bar-np.sqrt(nu_bar**(2) - (gamma_1*gamma_2*(1/18))))/ (0.5*gamma_2)
    collision_time_positive_bar= (nu_bar+np.sqrt(nu_bar**(2) - (gamma_1*gamma_2*(1/18))))/ (0.5*gamma_2)
    return collision_time_negative_bar,collision_time_positive_bar

collision_time_negative_bar=collision_time(nu_bar,gamma_1,gamma_2)[0]
collision_time_positive_bar=collision_time(nu_bar,gamma_1,gamma_2)[1]

md_timestep=collision_time_negative_bar/10 # just  change this to make the collision ratio integer
collision_ratio= collision_time_negative_bar/md_timestep
print("collision_ratio",collision_ratio)
    
box_size_bar= integer_box_size
print("box size to collison cell size ratio",box_size_bar/collision_cell_size_bar)

srd_count= (box_size_bar**3 )* fluid_particle_number_density

run_time_from_test= 26* 1.333# estimate taken from 4 procs and 2e6 time steps 
esimtated_run_time_mins= (srd_count/46305) *run_time_from_test
equilibration_timesteps_from_test=10000
estimated_equilibration_steps=(srd_count/46305) * equilibration_timesteps_from_test

nu= nu_bar * (lengthscale**2)/time_scale
Sc_est=62 # from pre calc 
diffusivity = nu_bar/Sc_est
characteristic_flow_vel=0.5
c_s = np.sqrt(5/3)
acoustic_timescale= box_size_bar/c_s
print(acoustic_timescale)
viscous_timescale= (box_size_bar**2)/nu_bar
print(viscous_timescale)
diffusion_timescale= (box_size_bar**2) / diffusivity # no embedded particles so doesnt matter
stokes_timescale= box_size_bar/characteristic_flow_vel

print("time_scale",time_scale)
print("lengthscale",lengthscale)
print("mass scale",mass_scale)
print("energy scale",energy_scale)
#print(collision_time_positive_bar)
#print(collision_time_negative_bar)
print("collision ratio",collision_ratio)
print("box size",box_size_bar)
print("srd count",srd_count)
print("estimated run time in hours", esimtated_run_time_mins/60)
print("estimated equilibration steps", estimated_equilibration_steps)

#%% Plots 

# plotting delta t vs nu bar 
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rcParams.update({'font.size': 25})
colour = [
 'black',
 'blueviolet',
 'cadetblue',
 'chartreuse',
 'coral',
 'cornflowerblue',
 'crimson',
 'darkblue',
 'darkcyan',
 'darkgoldenrod',
 'darkgray']
#nu_bar=np.linspace(np.sqrt(7/90),0.6,50)
#
# need to change this if you change fluid particle number density 
def gamma_1(M):
    gamma_1=1-((1-np.exp(-M))/M)
    return gamma_1
def gamma_2(M):
    gamma_2= (M+2)/(M-1)
    return gamma_2


M=np.array([5,100,1000])
gamma_1= gamma_1(M)
gamma_2= gamma_2(M)
gamma_1_2=gamma_1* gamma_2 *(1/18)

nu_bar= np.array(np.linspace(np.sqrt(gamma_1_2),1,50))

#plotting timestep vs nu_bar
def collision_time(nu_bar,gamma_1,gamma_2):
    collision_time_negative_bar= (nu_bar-np.sqrt(nu_bar**(2) - (gamma_1*gamma_2*(1/18))))/ (0.5*gamma_2)
    collision_time_positive_bar= (nu_bar+np.sqrt(nu_bar**(2) - (gamma_1*gamma_2*(1/18))))/ (0.5*gamma_2)
    return collision_time_negative_bar,collision_time_positive_bar

collision_time_negative_bar=collision_time(nu_bar,gamma_1,gamma_2)[0]
collision_time_positive_bar=collision_time(nu_bar,gamma_1,gamma_2)[1]

for i in range(M.size):
    
    plt.plot(nu_bar[:,i],collision_time_negative_bar[:,i],linestyle='dashed', label="$M(-)="+str(M[i])+"$", color=colour[i])
    
    plt.plot(nu_bar[:,i],collision_time_positive_bar[:,i],linestyle='dotted',label="$M(+)="+str(M[i])+"$", color=colour[i])
    plt.xlabel("$\\bar{\\nu}$")
    plt.ylabel("$\Delta \\bar{t}$",rotation=0, labelpad=10)
    plt.yscale('log')
    

plt.hlines(0.1,0.01 , 1, colors=None, linestyles='solid', label='$\mathrm{Kn}=0.1$',color='red')
plt.legend(loc='best', bbox_to_anchor=(1,1.05))
plt.savefig("scaled_timstep_vs_scaled_nu_M_5_.pdf",dpi=500, bbox_inches='tight')
plt.show()


#plotting viscosity contributions 
lamda_neg=collision_time_negative_bar
def kinetic_viscosity(lamda_neg,rotation_angle,M):
    angle_term= 1/(4-2*np.cos(rotation_angle)-2*np.cos(2*rotation_angle))
    M_term= (5*M)/(M-1)
    return lamda_neg * ((angle_term*M_term)-0.5)
def collisional_viscosity(lamda_neg,rotation_angle,M):
    angle_term=(1-np.cos(rotation_angle))/18
    M_term= (1-(1/M))
    return (1/lamda_neg)*angle_term*M_term
rotation_angle=np.pi/2
collisional_viscosity=collisional_viscosity(lamda_neg,rotation_angle,M)
kinetic_viscosity=kinetic_viscosity(lamda_neg,rotation_angle,M)
total_viscosity=collisional_viscosity+kinetic_viscosity
for i in range(1):
    
    plt.plot(lamda_neg[:,i],collisional_viscosity[:,i],linestyle='dashed',label="$\\bar{\\nu}_{coll},\  M="+str(M[i])+"$",color=colour[i])
    plt.plot(lamda_neg[:,i],kinetic_viscosity[:,i],linestyle='dotted', label="$\\bar{\\nu}_{kin},\ M="+str(M[i])+"$", color=colour[i])
   
  #  plt.plot(nu_bar[:,i],collision_time_positive_bar[:,i],label="$M(+)="+str(M[i])+"$")
    plt.xlabel("$\\bar{\lambda}$",labelpad=15)
    plt.ylabel("$ \\bar{\\nu}$",rotation=0, labelpad=20)
    plt.yscale('log')
    

plt.vlines(0.1,0 , 1, colors='red', linestyles='solid', label='$\mathrm{Kn}=0.1$')
plt.legend(loc='best', bbox_to_anchor=(1,1.05))
plt.savefig("predicted_eta_vs_scaled_mfp_.pdf",dpi=500, bbox_inches='tight')
plt.show()



 
#%% KATHLEEN paths 
##################
Path_2_shell_scirpts='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_KATHLEEN'
abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps-23Jun2022/build_kathleen/lmp_mpi'
abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/pure_srd_fix_deform.file'
Path_2_generic='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_KATHLEEN'
data_transfer_instructions=''
extra_code=''
wd_path='/home/ucahlrl/Scratch/output/'
num_task_req=''
tempdir_req=''
##################
### Variables settings
# calculating possible SRD_MD ratios and MD timesteps 
SRD_MD_ratio_ = 10
md_timestep=collision_time_negative_bar/SRD_MD_ratio_
#no_timestep_=(10000000*(SRD_MD_ratio_/SRD_MD_ratio_[0])).astype('int')
no_timestep_=np.array([10000])
print(md_timestep)
print(collision_time_negative_bar)
dump_freq_=10
thermo_freq=1000
erate_=np.array([0])
VP_ave_freq=10000
i_=0
j_=10
Path_2_generic='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_KATHLEEN'
hypthread='2'
fluid_name='mpcpure'
max_cores= 80
box_size_index=4


realisation_index_=[1,2,3,4,5,6,7,8,9,10]

ram_requirement='2G'
wall_time='2:00:00'
np_req=str(max_cores)
phi_= str(phi[box_size_index])
 # for one KATHLEEN node 
if (int(np_req)) > max_cores:
      print("Too many cores requested")
      breakpoint()
else:
      print("Core request satisfactory, producing simulation submission script ")

var_choice_1=erate_
var_choice_2=no_timestep_

def sim_file_prod_fix_deform_pure_individual_kathleen(nu_bar,erate_,var_choice_1,var_choice_2,hypthread,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name,SRD_MD_ratio_,dump_freq_):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_pure_fix_deform_mpc_visc_'+str(nu_bar)+'_box_size_'+str(box_size_bar[box_size_index])+'_gdot_'+str(erate_[0])+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    # test to check consistency of cores request 
    
    
    for j in range(i_,j_): #or now just use one realisation 
        for k in range(0,var_choice_1.size):       
                for m in range(0,var_choice_2.size):#range(0,1):  
                        param_set_code=str(np.random.randint(0, 1000000))
                        simulation_run_name=fluid_name+'_pure_fix_deform_mpc_visc_'+str(nu_bar)+'_'+str(sim_batchcode)+'_'+param_set_code+'_realisation_'+str(j)+'_erate_'+str(var_choice_1[k])+'_timestep_'+str(md_timestep)+'_no_timesteps_'+str(no_timestep_[m])+'_SRDMDratio_'+str(SRD_MD_ratio_)+'_'
                        run_code=''
                        no_SRD=str(int(srd_count[box_size_index])) 
                        #print(no_SRD)
                        box_size = str(box_size_bar[box_size_index])
                        timestep_input= str(md_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=str(SRD_MD_ratio_)
                        lamda= str(collision_time_negative_bar)
                        dump_freq=str(dump_freq_) # 
                        thermo_freq = str(thermo_freq) # 
                        temp_=str(1)
                        grid_size=str(1)
                        mass_SRD=str(1)
                        no_timesteps = str(no_timestep_[m])
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))
                        num_proc=str(np_req)
                        erate=str(erate_[m])
                        


                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var erate_in '+str(erate)+'  -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+'  -var VP_ave_freq '+str(VP_ave_freq)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+'  -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -in '+abs_path_2_lammps_script+' \n &'  #>> '+prod_run_file_name+' & \n'
                                            # -var erate_in 0 -var fluid_name equilibrium -var  sim_batchcode 806324 -var swap_rate 15 -var swap_number 1 -var VP_ave_freq 10000 -var equilibration_timesteps 1000 -var chunk 20 -var grid_size 1 -var realisation_index 2 -var temp_ 1 -var lambda 0.05071624521210362 -var rand_int 1234 -var rand_int_1 5678 -var rand_int_2 13565 -var no_SRD 60835 -var mass_SRD 1 -var box_size 23.0 -var timestep_input 0.005071624521210362 -var SRD_MD_ratio 10 -var dump_freq 1000 -var thermo_freq 1000 -var no_timesteps 600000 -var vel_target 10


                        #print(run_code_individual)
                        run_code=run_code +run_code_individual


                        run_code = run_code[:-3]

                        py2bash_launch_overwriter_kathleen.py2bash_launch_overwriter(hypthread,Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)


# %%
sim_file_prod_fix_deform_pure_individual_kathleen(nu_bar,erate_,var_choice_1,var_choice_2,hypthread,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name,SRD_MD_ratio_,dump_freq_)
       
# %%
