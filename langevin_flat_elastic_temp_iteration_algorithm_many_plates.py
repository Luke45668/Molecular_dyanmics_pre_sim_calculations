#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 5/12/23

This script will does all the pre calcs for the SRD fluid mappings and flat elastic particles then produces the simulation run files for Kathleen. 


@author: lukedebono
"""
#%%
#from msilib import MSIMODIFY_INSERT_TEMPORARY
import os
import numpy as np
import sigfig as sgf
import matplotlib.pyplot as plt
import regex as re
import pandas as pd
import math as m 
import sigfig 
from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
import seaborn as sns
import math as m
import scipy.stats
from datetime import datetime
from sim_file_producer_SRD import *
import subprocess as subproc
import sys



#%% Producing random dist of orientations 
# produce N random numbers between 0,1 
box_size_bar=np.array([100])

number_of_points=5 # for bug test , normally 30
# could put in spencers PRNG 







# %% producing shell scripts for MYRIAD with rotations 
# on my computer 


Path_2_shell_scirpts='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'

# for running on my computer 
abs_path_2_lammps_exec='/Users/luke_dev/Documents/lammps_hirotori/build_m3maxbook/lmp_mpi'

abs_path_2_lammps_script='/Users/luke_dev/Documents/LSC/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_mol_rattle' 


Path_2_generic='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'


num_task_req=''
data_transfer_instructions=''
SRD_MD_ratio_ = 10
VP_ave_freq=10000
md_timestep=0.005071624521210362
collision_time_negative_bar=0.05071624521210362


#erate= np.array([0.0075,0.005,0.0025]) longer runs which need checkpointing
erate=np.array([1,0.9,0.7,0.5,0.35,0.2,0.1,0.09,0.08,
                0.07,0.06,0.05,0.04,
                0.03,0.0275,0.025,0.0225,
                0.02,0.0175,0.015,0.0125,
                0.01,0.0075,0.005,0.0025,
                0.001,0.00075,0.0005])
# erate=np.array([100,50,25,10,5])

i_=0
j_=number_of_points
fluid_name='langevinrun'
bending_stiffness=np.array([1000]) # original 50,100,200,400


damp_init=0.03633 # based on H20 as reference fluid 


internal_stiffness=np.array([500])





#no_timestep_=np.array([2000000]) 
# just for dist test 
#no_timestep_=np.array([1000]) 
#phantom_mass=np.array([0.01])
equilibrium_triangle_side_length=3




realisation_index_=np.arange(0,1000,1)
timestep_multiplier=0.2

def compute_timesteps_for_strain(total_strain,erate,md_timestep,timestep_multiplier):
      no_timestep_=(np.round((total_strain/(erate*md_timestep*timestep_multiplier)),decimals=-3)).astype('int')

      return no_timestep_

# for MYRIAD run
total_strain=400
no_timestep_=compute_timesteps_for_strain(total_strain,erate,md_timestep,timestep_multiplier)
if np.any(no_timestep_>2e9):
     print("error! too many timesteps, must be less than 2e9")



# # for bug test
# number_of_restarts_per_run=np.array([1,1,1,1,1])
# no_timestep_=np.array([2000,2000,2000,2000,2000])
# np_req=np.array([14,14,14,14,14]).astype('str')

def folder_check_or_create(filepath,folder):
     os.chdir(filepath)
     # combine file name with wd path
     check_path=filepath+"/"+folder
     print((check_path))
     if os.path.exists(check_path) == 1:
          print("file exists, proceed")
          os.chdir(check_path)
     else:
          print("file does not exist, making new directory")
          os.chdir(filepath)
          os.mkdir(folder)
          os.chdir(filepath+"/"+folder)


# initial coordinates 
j=0
n=0
l=0 # intstiff parameter
# data so far 

input_temp_seed =np.array([
 0.75,
 0.75,
 0.85,
 0.9,
 0.925,
 0.95,
1,
 1,
 1,
 1,
 1,
1,
 1,
 1,
 1,
 1,
1,
 1,
 1,
 1,
 1,
1,
1, 
1,
1,
1,
1,
1])
count=0
tolerance=0.05
increment=0.01
input_temp=1
desired_temp=1
from log2numpy import *
thermo_vars='         KinEng         PotEng         Press         c_myTemp        c_bias         TotEng    '
filepath='/Users/luke_dev/Documents/simulation_run_folder/many_plates_production_run_rattle'
os.chdir(filepath)

new_temp=[]
new_temp_array=np.zeros((erate.size))
dump_freq=str(100)
thermo_freq=str(100)
temp_comparison_list_conv=[]
input_temp_list_conv=[]
#%%
# this could be parallelised probably qutie easily

# need to use this loop to find a temp which gives T=1
# then another loop to run it until we get 10 sucessful runs 
erate_1=int(sys.argv[1]) 
erate_2=int(sys.argv[2])
for n in range(erate_1,erate_2):
    input_temp=input_temp_seed[n]
    

    temp_comparison=0.5
    for j in range(1):
    
        #input_temp=input_temp_seed[n]
        temp_comparison=0.5
        while np.abs(temp_comparison)>tolerance:
        # for 8,16,32,36
            var_choice_1=erate
            var_choice_2=internal_stiffness
            # individual shear rates 

            run_code_list=[]
            # test to check consistency of cores request 


            run_code=''
            erate_in=erate[n]
            #print(no_SRD)
            box_size = str(100)
            timestep_input= str(md_timestep)
            # number of chunks to use for VP averaging
            SRD_MD_ratio=str(int(SRD_MD_ratio_))
            lamda= str(collision_time_negative_bar)

            no_timesteps = str(no_timestep_[n])
            rand_int =str(np.random.randint(0, 1000000))
            rand_int_1 =str( np.random.randint(0, 1000000))
            rand_int_2 =str(np.random.randint(0, 1000000))
            rand_int_3=str(np.random.randint(0,1000000))

                        
            run_code_individual ="mpirun -np 2 "+abs_path_2_lammps_exec+' -var temp '+str(input_temp)+' -var damp '\
                +str(damp_init)+' -var erate_in '+str(erate_in)+' -var equilirbium_triangle_side_length '\
                    +str(equilibrium_triangle_side_length)+\
            ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[l])+\
                ' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(999)+\
            ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[j])+\
                ' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
            ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var box_size '+box_size+\
                ' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
            ' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+\
                no_timesteps+' -in '+abs_path_2_lammps_script+' \n '  #>> '+prod_run_file_name+' & \n'

            #os.system(run_code_individual)
            subproc.run(run_code_individual,shell=True)

            logfile_name="log.langevinrun_no"+str(999)+"_hookean_flat_elastic_"+rand_int+\
                "_"+str(realisation_index_[j])+"_"+str(box_size)+"_0.03633_0.005071624521210362_"+dump_freq+"_"+thermo_freq+"_"+no_timesteps+"_0.2_gdot_"+str(erate_in)+\
                "_BK_1000_K_"+str(internal_stiffness[l])
            try: 
                log_file_data=log2numpy_reader(logfile_name,
                                        filepath,
                                        thermo_vars)
            except:
                continue
            
            

            temp_out=np.mean(log_file_data[50:,5])
            temp_comparison=desired_temp-temp_out
            temp_comparison_list_conv.append(temp_comparison)
            input_temp_list_conv.append(input_temp)
            print("Temp comparison", temp_comparison)
            count+=1



            if np.abs(temp_comparison)>tolerance:
                if temp_comparison>0: # pos
                    input_temp+=increment

                    if input_temp<0:
                     breakpoint
                       
                     print("New input temp",input_temp)
                elif temp_comparison<0:
                    input_temp-=increment
                    if input_temp<0:
                     breakpoint
                    
                     print("New input temp",input_temp)

            else:
                print("temp controlled")
                new_temp_array[n]=input_temp
                os.system("cp -r "+logfile_name+ " passed_tests/") #not working 
                os.system("cp -r *_"+rand_int+"_"+str(realisation_index_[j])+"*h5 passed_tests/")
                os.system("cp -r *_"+rand_int+"_"+str(realisation_index_[j])+"*dump passed_tests/")


#%% run at fixed t

for n in range(erate_1,erate_2):
    count=1
    
    
    while count < 5:
    
        #input_temp=input_temp_seed[n]
        temp_comparison=0.5
        while np.abs(temp_comparison)>tolerance:
        # for 8,16,32,36
            var_choice_1=erate
            var_choice_2=internal_stiffness
            # individual shear rates 

            run_code_list=[]
            # test to check consistency of cores request 


            run_code=''
            erate_in=erate[n]
            #print(no_SRD)
            box_size = str(100)
            timestep_input= str(md_timestep)
            # number of chunks to use for VP averaging
            SRD_MD_ratio=str(int(SRD_MD_ratio_))
            lamda= str(collision_time_negative_bar)

            no_timesteps = str(no_timestep_[n])
            rand_int =str(np.random.randint(0, 1000000))
            rand_int_1 =str( np.random.randint(0, 1000000))
            rand_int_2 =str(np.random.randint(0, 1000000))
            rand_int_3=str(np.random.randint(0,1000000))

                        
            run_code_individual ="mpirun -np 2 "+abs_path_2_lammps_exec+' -var temp '+str(input_temp)+' -var damp '\
                +str(damp_init)+' -var erate_in '+str(erate_in)+' -var equilirbium_triangle_side_length '\
                    +str(equilibrium_triangle_side_length)+\
            ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[l])+\
                ' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(999)+\
            ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[int(count)])+\
                ' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
            ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var box_size '+box_size+\
                ' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
            ' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+\
                no_timesteps+' -in '+abs_path_2_lammps_script+' \n '  #>> '+prod_run_file_name+' & \n'

            #os.system(run_code_individual)
            subproc.run(run_code_individual,shell=True)

            logfile_name="log.langevinrun_no"+str(999)+"_hookean_flat_elastic_"+rand_int+\
                "_"+str(realisation_index_[int(count)])+"_"+str(box_size)+"_0.03633_0.005071624521210362_"+dump_freq+"_"+thermo_freq+"_"+no_timesteps+"_0.2_gdot_"+str(erate_in)+\
                "_BK_1000_K_"+str(internal_stiffness[l])
            try: 
                log_file_data=log2numpy_reader(logfile_name,
                                        filepath,
                                        thermo_vars)
            except:
                continue
            

            temp_out=np.mean(log_file_data[5:,5])
            temp_comparison=desired_temp-temp_out
            temp_comparison_list_conv.append(temp_comparison)
            input_temp_list_conv.append(input_temp)
            print("Temp comparison", temp_comparison)
            

         
               

            if np.abs(temp_comparison)<tolerance:
                print("temp controlled")
                count+=1
                
                os.system("cp -r "+logfile_name+ " passed_tests/") #not working 
                os.system("cp -r *"+rand_int+"*h5 passed_tests/")
                os.system("cp -r *"+rand_int+"*dump passed_tests/")







     
     
        


# %%
plt.plot(temp_comparison_list_conv,label="temp comp")
plt.plot(input_temp_list_conv,label="input T")
plt.axhline(0.05)
plt.axhline(-0.05)
plt.legend()
plt.show()
# %%
