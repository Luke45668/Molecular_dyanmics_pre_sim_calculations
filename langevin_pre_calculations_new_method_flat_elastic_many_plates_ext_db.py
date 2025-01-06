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
from simulation_production_module import *



box_size_bar=100


number_of_points=20



Path_2_shell_scirpts='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'

# for running on my computer 
abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps_hirotori/build_MYRIAD_ext/lmp_mpi'
#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_mol'

#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/lammps_scripts/in.nvt_brownian_uef_flat_elastic_particles'
#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/lammps_scripts/in.nvt_brownian_uef_flat_elastic_particles_biax'
abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/lammps_scripts/in.nvt_uef_oldroyd_db_ext'


#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/lammps_scripts/in.nvt_uef_oldroyd_db'

#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/lammps_scripts/in.brownian_uef_flat_elastic_particles_biax'

#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.langevin_with_hookean_flat_elastic_particle_pentagon_mol_rattle'
#for running on myriad 
# abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps_hirotori/build_serial/lmp'
# abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_eq_b4_shear'

Path_2_generic='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'
extra_code='module unload mpi compilers gcc-libs \n module load beta-modules \n module load gcc-libs/10.2.0 \n module load compilers/intel/2022.2 \n module load mpi/intel/2019/update6/intel \n  module load hdf/5-1.12.3-impi/intel-2022'
wd_path='/home/ucahlrl/Scratch/output/nvt_runs/final_plate_runs/'


wd_path='/home/ucahlrl/Scratch/output/nvt_runs/db_run_100_small_tstep_novisc'
num_task_req=''
data_transfer_instructions=''
SRD_MD_ratio_ = 10
VP_ave_freq=10000
md_timestep=0.005071624521210362
collision_time_negative_bar=0.05071624521210362


erate=np.linspace(1,0.005,24)
erate=np.logspace(-5,-0.25,24)

i_=0
j_=number_of_points

fluid_name='DBextnvt'


bending_stiffness=np.array([500]) 

internal_stiffness=np.array([100,150,300,600])
internal_stiffness=np.array([30,60,100,150,300,600])
internal_stiffness=np.array([5,7.5,10,15,20,30])
internal_stiffness=np.array([5,15,30,60,90,120])
#internal_stiffness=np.array([100,200])

# internal_stiffness=np.array([30,60])

# internal_stiffness=np.array([120,240])

#internal_stiffness=np.array([480,960])

#internal_stiffness=np.array([1250,1500])

#internal_stiffness=np.array([3000,6000])


damp=np.array([0.035])


equilibrium_triangle_side_length=3
tempdir_req='1G'
ram_requirement='5G'
wall_time='48:00:00'


# input_temp=np.array([0.75,0.80,
#  0.85,0.9,0.925,0.95,0.95,0.95,0.95,0.95,
# 0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975,0.975])

# input_temp=np.array([0.8,0.825,0.85,
#  0.875,0.9,0.925,0.95,0.975,1,1,1,
# 1,1,1,1,1,1,1,1,1])

# no langevin temp
input_temp=np.array([1,1,1,
 1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])


realisation_index_=np.arange(0,1000,1)
# need to make this an array for each
#timestep_multiplier=0.2

# for flat elastic 
timestep_multiplier=np.array([
[0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.0005,0.0005,0.0005,
0.0005,0.0005,0.005,0.005,0.2],

[0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.0005,0.0005,0.0005,
0.0005,0.0005,0.005,0.005,0.2]])

# for dumbell test 
timestep_multiplier=np.array([
[0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,
0.0005,0.0005,0.0005,0.0005,0.0005,0.005,
0.005],

[0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,
0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,
0.0005,0.0005,0.0005,0.0005,0.0005,0.005,
0.005]])*4

# for plate test 
timestep_multiplier=np.array([[3.94351157e-05, 4.12088897e-05, 4.31497461e-05, 4.52824594e-05,
       4.76369570e-05, 5.02497320e-05, 5.31657481e-05, 5.64410493e-05,
       6.01463967e-05, 6.43724388e-05, 6.92372261e-05, 7.48974122e-05,
       8.15654372e-05, 8.95367880e-05, 9.92349740e-05, 1.11289284e-04,
       1.26677048e-04,0.002 , 0.002 , 0.002 , 0.002 , 0.002 , 0.02  , 0.02  ],
       [0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.002 , 0.002 , 0.002 , 0.002 , 0.002 , 0.02  , 0.02  ],
       [0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.002 , 0.002 , 0.002 , 0.002 , 0.002 , 0.02  , 0.02  ],
       [0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.002 , 0.002 , 0.002 , 0.002 , 0.002 , 0.02  , 0.02  ],
       [0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002,
        0.0002, 0.002 , 0.002 , 0.002 , 0.002 , 0.002 , 0.02  , 0.02  ]])

thermal_damp_multiplier=np.array([25,25,25,25,25,25,25,100,100,100,100,100,
100,100,100,100,250,250])/10

# thermal_damp_multiplier=np.array([100,100,100,100,100,100,100,100,100,100,100,100,
# 100,100,100,100,150,150])/10
total_strain=100
def min_timestep_multi(total_strain,erate,md_timestep):
    
    min_multi=total_strain/(erate*(2e9-1000)*md_timestep)

    return min_multi

min_multi=min_timestep_multi(total_strain,erate,md_timestep)
# for MYRIAD run
min_multi=np.tile(min_multi,(internal_stiffness.size,1))
timestep_multiplier=min_multi
no_timestep_=compute_timesteps_for_strain(total_strain,erate,md_timestep,min_multi)
#no_timestep_[:,-1]=10000000
if np.any(no_timestep_>2e9):
     print("error! too many timesteps, must be less than 2e9")
#no_timestep_[:,-1]=10000000 #make equilibrium 10 mil steps 

# only for equilibrium run 

# dump_freq=out_put_freq_calc(no_timestep_,10000)
# thermo_freq=out_put_freq_calc(no_timestep_,10000)

dump_freq=out_put_freq_calc(no_timestep_,1000)
thermo_freq=out_put_freq_calc(no_timestep_,1000)


np_req=str(8) 
var_choice_1=erate
var_choice_2=internal_stiffness

thermal_damp_multiplier=np.repeat(250,24)
thermal_damp_time=timestep_multiplier*md_timestep*thermal_damp_multiplier
timestep=timestep_multiplier*md_timestep

#NOTE: for next set of runs change in.nvt_uef so that stretch is along z axis 
#NOTE: comms cut off is only 6, could push it up to 10 at some cost 

# %%
# all in one file for myriad
# run tests with various sets of t damp
    
sim_file_prod_flat_elastic_MYRIAD_all_erate_one_file(SRD_MD_ratio_,collision_time_negative_bar,
                                                    bending_stiffness,damp,
                                                    input_temp,
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
                                                    fluid_name,timestep_multiplier,thermal_damp_multiplier)
        


#%% fix a single t damp to each shear rate

sim_file_prod_flat_elastic_MYRIAD_all_erate_one_file_var_damp(SRD_MD_ratio_,collision_time_negative_bar,bending_stiffness,damp,input_temp,
                                                         erate,
                                                         equilibrium_triangle_side_length,
                                                         var_choice_1,var_choice_2,internal_stiffness,
                                                         data_transfer_instructions,extra_code,wd_path,
                                                         np_req,num_task_req,tempdir_req,wall_time,
                                                         ram_requirement,realisation_index_,VP_ave_freq,
                                                         abs_path_2_lammps_exec,abs_path_2_lammps_script,
                                                         no_timestep_,thermo_freq,dump_freq,md_timestep,
                                                         i_,j_,box_size_bar,Path_2_shell_scirpts,Path_2_generic,fluid_name,timestep_multiplier,thermal_damp_multiplier)


# %% producing files in simulation test folder only for equilibrium 
abs_path_2_lammps_exec="/Users/luke_dev/Documents/lammps_hirotori/build_m3maxbook/lmp_mpi"

abs_path_2_lammps_script="/Users/luke_dev/Documents/LSC/langevin_codes/shear_bending_stiff_approx_codes/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_mol_rattle"

os.chdir("/Users/luke_dev/Documents/simulation_run_folder/")
h=2 # damp=0.035
path="/Users/luke_dev/Documents/simulation_run_folder/eq_run_tri_plate_damp_"+str(damp[h])+"_K_500_4000"
# os.mkdir(path)
os.chdir(path)
run_code_list=[]
h= # damp=0.035
l=0 #K=4000, damp=0.035
input_temp=np.array([1])
erate=np.array([0])
no_timestep_=np.array([1000000])
os.chdir(path)
MyFile=open('K_'+str(internal_stiffness[l])+'_damp_'+str(damp[h])+'_erate_'+str(erate[0])+'.sh','w')
sim_batchcode=str(np.random.randint(0, 1000000))


for k in range(1):    
            
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
           
for element in run_code_list:
            MyFile.write(element)
            MyFile.write('\n')
MyFile.close()
# %%
