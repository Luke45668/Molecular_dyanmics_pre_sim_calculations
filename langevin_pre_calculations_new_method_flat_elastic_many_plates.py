#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 5/12/23

This script will does all the pre calcs flat elastic particle simulations from a template then produces the simulation run files for Kathleen. 


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



box_size_bar=100


number_of_points=20







# %% producing shell scripts for MYRIAD with rotations 
# on my computer 
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_imac_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_macbook_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_macbook_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_m3maxbook/lmp_mpi'
#abs_path_2_lammps_script='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/LSC/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5'
#abs_path_2_lammps_script='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/LSC/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5_chain'


Path_2_shell_scirpts='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'

# for running on my computer 
abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps_hirotori/build_MYRIAD/lmp_mpi'
abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_mol'

#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_mol_rattle'

#for running on myriad 
# abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps_hirotori/build_serial/lmp'
# abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_eq_b4_shear'


Path_2_generic='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'



extra_code='module unload mpi compilers gcc-libs \n module load beta-modules \n module load gcc-libs/10.2.0 \n module load compilers/intel/2022.2 \n module load mpi/intel/2019/update6/intel \n  module load hdf/5-1.12.3-impi/intel-2022'
wd_path='/home/ucahlrl/Scratch/output/'
num_task_req=''
data_transfer_instructions=''
SRD_MD_ratio_ = 10
VP_ave_freq=10000
md_timestep=0.005071624521210362
collision_time_negative_bar=0.05071624521210362
#erate= np.array([0.0008,0.001,0.002,0.005,0.01,0.1])
erate= np.array([0.03,0.0275,0.025,0.0225,0.02,0.0175,0.015,0.0125,0.01,0.0075])
#erate=np.array([0.005,0.0025,0.001,0.00075,0.0005]) 


#erate= np.array([0.0075,0.005,0.0025]) longer runs which need checkpointing
erate=np.array([1,0.9,0.7,0.5,0.2,0.1,0.09,0.08,
                0.07,0.06,0.05,0.04,
                0.03,0.0275,0.025,0.0225,
                0.02,0.0175,0.015,0.0125,
                0.01,0.0075,0.005,0.0025,
                0.001,0.00075,0.0005])

# erate=np.array([1,0.8,0.6,0.4,0.2,0.1,0.08,
#                 0.06,0.04,
#                 0.03,0.025,
#                 0.02,0.015,
#                 0.01,0.005,
#                 0.001,0.0005])
# erate=np.array([100,50,25,10,5])

i_=0
j_=number_of_points
fluid_name='langevinrun'
bending_stiffness=np.array([500]) # original 50,100,200,400
bending_stiffness=np.array([10000])
#internal_stiffness=np.array([60,80,100]) # 20 does nothing 
def product_constraint_inputs(damp_upper_bound,damp_lower_bound,internal_stiff,initial_damp,n_sets):
      kdamp=internal_stiff*initial_damp
      Possible_damps=np.linspace(damp_lower_bound,damp_upper_bound,n_sets)

      new_internal_stiffness=kdamp/Possible_damps

      return new_internal_stiffness,Possible_damps

      # fix 


# damp_init=0.03633 # based on H20 as reference fluid 
# damp_init=0.05
# new set  with fixed product=kdamp
#internal_stiffness,damp=product_constraint_inputs(0.1,0.01,internal_stiffness_init,damp_init,10)


internal_stiffness=np.array([500,2000])

damp=np.array([0.01,0.035,0.05])




#no_timestep_=np.array([2000000]) 
# just for dist test 
#no_timestep_=np.array([1000]) 
#phantom_mass=np.array([0.01])
equilibrium_triangle_side_length=3
tempdir_req='1G'
ram_requirement='10G'
wall_time='48:00:00'


np_req=np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]).astype('str')

# inital_temp_multipl=np.array([1.06006728, 0.99582948, 0.98558312, 1.01943334, 1.01784   ,
#        0.94875864, 0.9455115 , 0.88103351, 0.9308963 , 0.86103214])



input_temp=np.array([0.75,0.75,
 0.85,0.9,0.925,0.95,0.95,0.95,0.95,0.95,
0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975, 
0.975,0.975,0.975,0.975,0.975])

realisation_index_=[1,2,3]
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





dump_freq=np.array([100, 100, 100, 100, 100, 100,1000, 1000,
 1000, 1000, 1000, 1000, 1000, 1000, 1000,
 1000, 10000, 10000, 10000, 10000, 10000, 10000,
 10000, 10000, 10000, 10000, 10000])
thermo_freq=np.array([100, 100, 100, 100, 100, 100,1000, 1000,
 1000, 1000, 1000, 1000, 1000, 1000, 1000,
 1000, 10000, 10000, 10000, 10000, 10000, 10000,
 10000, 10000, 10000, 10000, 10000])
np_req=str(8)
var_choice_1=erate
var_choice_2=internal_stiffness
# individual shear rates 


# %%
# all in one file for myriad
def sim_file_prod_flat_elastic_MYRIAD_all_erate_one_file(damp,input_temp,
                                                         erate,
                                                         equilibrium_triangle_side_length,
                                                         var_choice_1,var_choice_2,internal_stiffness,
                                                         data_transfer_instructions,extra_code,wd_path,
                                                         np_req,num_task_req,tempdir_req,wall_time,
                                                         ram_requirement,realisation_index_,VP_ave_freq,
                                                         abs_path_2_lammps_exec,abs_path_2_lammps_script,
                                                         no_timestep_,thermo_freq,dump_freq,md_timestep,
                                                         i_,j_,box_size_bar,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_eqts_'+\
        str(equilibrium_triangle_side_length)+'_realisations_'+str(j_)+'_box_size_'+\
            str(box_size_bar)+'_bendstiff_'+str(bending_stiffness[0])+\
                '_intstiff_'+str(internal_stiffness[0])+'_'+str(internal_stiffness[-1])+\
                    '_erate_'+str(var_choice_1[0])+'_'+str(var_choice_1[-1])+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    run_code_list=[]
    # test to check consistency of cores request 

    

    #for n in range(0,np_req.size):
        #or now just use one realisation 
    for h in range(damp.size):

        for k in range(erate.size):    
        #for k in range(0,1):   
            for m in range(0,var_choice_2.size): 
                        #for m in range(0,1):  
                                for j in range(i_,j_):
                                    param_set_code=str(np.random.randint(0, 1000000))
                                    simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+\
                                        '_realisation_'+str(j)+'_Bk_'+str(bending_stiffness[0])+'_np_'+str(np_req)+\
                                            '_no_timesteps_'+str(no_timestep_[k])+'_intstiff_'+str(var_choice_2[m])+\
                                                '_eqsl_'+str(equilibrium_triangle_side_length)+'_erate_'+\
                                                    str(var_choice_1[k])+'_damp_'+str(damp[h])+'_'
                                    run_code=''
                                
                                    #print(no_SRD)
                                    box_size = str(box_size_bar)
                                    timestep_input= str(md_timestep)
                                    SRD_MD_ratio=str(int(SRD_MD_ratio_))
                                    lamda= str(collision_time_negative_bar)
                                    dump_freq_=str(dump_freq[k])
                                    thermo_freq_ = str(thermo_freq[k])
                                    no_timesteps = str(no_timestep_[k])
                                    rand_int =str(np.random.randint(0, 1000000))
                                    rand_int_1 =str( np.random.randint(0, 1000000))
                                    rand_int_2 =str(np.random.randint(0, 1000000))
                                    rand_int_3=str(np.random.randint(0,1000000))
                                    
                                

                                                
                                    run_code_individual ="mpirun -np "+np_req+" "+abs_path_2_lammps_exec+' -var temp '+str(input_temp[k])+' -var damp '\
                                        +str(damp[h])+' -var erate_in '+str(erate[k])+' -var equilirbium_triangle_side_length '\
                                            +str(equilibrium_triangle_side_length)+\
                                    ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[m])+\
                                        ' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+\
                                    ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[j])+\
                                        ' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
                                    ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var box_size '+box_size+\
                                        ' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
                                    ' -var dump_freq '+dump_freq_+' -var thermo_freq '+thermo_freq_+' -var no_timesteps '+\
                                        no_timesteps+' -in '+abs_path_2_lammps_script+' \n '  #>> '+prod_run_file_name+' & \n'

                                    
                                    run_code_list.append(run_code_individual)
                                    run_code=run_code +run_code_individual

                                    
                                    run_code = run_code[:-2]

                                    py2bash_launch_overwriter.py2bash_launch_overwriter_mpi(Path_2_generic,
                                                                                            simulation_batch_folder,
                                                                                            simulation_run_name,
                                                                                            specific_email,wall_time,
                                                                                            ram_requirement,tempdir_req,
                                                                                            num_task_req,np_req,wd_path,
                                                                                            extra_code,run_code,data_transfer_instructions)
    return run_code_list, sim_batchcode
    



sim_file_prod_flat_elastic_MYRIAD_all_erate_one_file(damp,input_temp,
                                                    erate,
                                                    equilibrium_triangle_side_length,
                                                    var_choice_1,
                                                    var_choice_2,
                                                    internal_stiffness,
                                                    data_transfer_instructions,
                                                    extra_code,wd_path,
                                                    np_req,
                                                    num_task_req,
                                                    tempdir_req,
                                                    wall_time,
                                                    ram_requirement,
                                                    realisation_index_,
                                                    VP_ave_freq,
                                                    abs_path_2_lammps_exec,
                                                    abs_path_2_lammps_script,
                                                    no_timestep_,
                                                    thermo_freq,
                                                    dump_freq,
                                                    md_timestep,
                                                    i_,
                                                    j_,
                                                    box_size_bar,
                                                    Path_2_shell_scirpts,
                                                    Path_2_generic,
                                                    fluid_name)
        




# %%
