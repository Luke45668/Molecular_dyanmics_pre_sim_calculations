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




#%% Producing random dist of orientations 
# produce N random numbers between 0,1 
box_size_bar=np.array([23])
#NOTE using mathematics convention for the angles
# r 
# theta is the equator 
#phi is the incline


number_of_points=100 # for bug test , normally 30
# could put in spencers PRNG 

def producing_random_points_with_theta(number_of_points):

    rng = np.random.default_rng(12345)
    Phi=np.arccos(1-2*(rng.random((number_of_points))))
    
    Theta=2*np.pi*rng.random((number_of_points))
    rho=1#7.7942286341
    A=Phi
    B=Theta
    R=np.array([rho*np.sin(A)*np.cos(B),rho*np.sin(B)*np.sin(A),rho*np.cos(A)])


    return Phi,Theta,R

producing_random_points_with_theta_out=producing_random_points_with_theta(number_of_points)
Phi=producing_random_points_with_theta_out[0]
Theta=producing_random_points_with_theta_out[1]
Rotation_vector=producing_random_points_with_theta_out[2]

plt.scatter(Theta,Phi)
pi_theta_ticks=[ 0, np.pi/2,np.pi, 3*np.pi/2, 2*np.pi]
pi_theta_tick_labels=['0', 'π/2', 'π', '3π/2', '2π']
plt.xticks(pi_theta_ticks, pi_theta_tick_labels)
plt.xlabel('$\\theta$')
pi_phi_ticks=[ 0,np.pi/4, np.pi/2,3*np.pi/4,np.pi]
pi_phi_tick_labels=[ '0','π/4', 'π/2','3π/4' ,'π']
plt.yticks(pi_phi_ticks,pi_phi_tick_labels)
plt.ylabel('$\phi$')

plt.show()

#checking on sphere plot
r=1
x=r*np.cos(Phi)*np.sin(Theta)
y=r*np.sin(Phi)*np.sin(Theta)
z=r*np.cos(Theta)
fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.scatter(x,y,z)
plt.show()

bin_count = int(np.ceil(np.log2((number_of_points))) + 1)# sturges rule
#bin_count=int(number_of_points/10 ) 
pi_theta_ticks=[ 0, np.pi/2,np.pi, 3*np.pi/2, 2*np.pi]
pi_theta_tick_labels=['0', 'π/2', 'π', '3π/2', '2π']
plt.hist(Theta, bins=bin_count,density=True)
plt.xlabel('$\\theta$')

plt.xticks(pi_theta_ticks, pi_theta_tick_labels)

plt.tight_layout()
#plt.savefig("theta_pre_histogram_"+str(number_of_points)+"_points.pdf",dpi=1200)

plt.show()

pi_phi_ticks=[ 0,np.pi/4, np.pi/2,3*np.pi/4,np.pi]
pi_phi_tick_labels=[ '0','π/4', 'π/2','3π/4' ,'π']
plt.hist(Phi, bins=bin_count,density=True)
plt.xlabel('$\phi$')
plt.xticks(pi_phi_ticks,pi_phi_tick_labels)
plt.tight_layout()
#plt.savefig("phi_pre_histogram_"+str(number_of_points)+"_points.pdf",dpi=1200)
plt.show()
#%% generating random orientations 

equilibrium_triangle_side_length=3


from mol_file_overwriter import *
from flat_elastic_vector_funcs_module import *

box_size_index=0

centre_of_plane=point_on_plane(
        box_size_bar,
        box_size_index,
                )



coordinates_tuple_3d=()
#j=50
for j in range(number_of_points):
    
    

    coordinates_2d=bead_general_coordinates_eq_triangle(
                                                            centre_of_plane[0],
                                                            centre_of_plane[1],
                                                            equilibrium_triangle_side_length
                                                                )
    
    #solve for bottom left stokes bead 
    z=solve_plane_equation_for_z(
            Rotation_vector[:,j],
            centre_of_plane,
            3,
            j, 
            coordinates_2d,
                )

    #print(z)

    phantom_bead_1=np.array([coordinates_2d[3][0],coordinates_2d[3][1],z])
   
    # calculate  basis vectors

    basis_vector_1=(phantom_bead_1-centre_of_plane)/compute_vector_magnitude(phantom_bead_1-centre_of_plane)   
    
    basis_vector_2= np.cross(Rotation_vector[:,j],(phantom_bead_1-centre_of_plane))/compute_vector_magnitude(np.cross(Rotation_vector[:,j],(phantom_bead_1-centre_of_plane)))
    
    # now compute new coordinates in terms of basis 

    coordinates_3d=equilateral_coords_from_basis(centre_of_plane, basis_vector_1,basis_vector_2,equilibrium_triangle_side_length)

    
    
    ell_1=coordinates_3d[1]-coordinates_3d[0]
    ell_1_mag=compute_vector_magnitude(ell_1)
    #print("ell_1_mag",ell_1_mag)
    ell_2=coordinates_3d[2]-coordinates_3d[0]
    ell_2_mag=compute_vector_magnitude(ell_2)
    #print("ell_2_mag",ell_2_mag)
    ell_3=coordinates_3d[2]-coordinates_3d[1]
    ell_3_mag=compute_vector_magnitude(ell_3)
    #print("ell_3_mag",ell_3_mag)
    ell_sum=ell_1_mag+ell_2_mag+ell_3_mag

    if m.isclose(9,ell_sum):
       #print(coordinates_3d)

       coordinates_tuple_3d=coordinates_tuple_3d+(coordinates_3d,)
    else: 
          print("calculation error")
          break 
          
    print("ell_1 cross ell_2",compute_vector_magnitude(np.cross(ell_1,ell_2)))

    # fig = plt.figure()
    # ax = plt.axes(projection ='3d')

    # for i in range(len(coordinates_2d)):
    #     x=coordinates_3d[i][0]
    #     y=coordinates_3d[i][1]
    #     z=coordinates_3d[i][2]
    #     ax.scatter(x,y,z)
    # plt.show()
    # print(coordinates_3d)

    # coordinates_tuple_3d=coordinates_tuple_3d+(coordinates_3d,)

    # now lets check the dist of ell_1 cross ell_2

#%% now lets re-check the distribution  of theta and phi
spherical_coordinates_area_vector=np.zeros((number_of_points,3))
#bin_count = int(np.ceil(np.log2((number_of_points))) + 1)
    


for i in range(number_of_points):
    ell_1=coordinates_tuple_3d[i][1]-coordinates_tuple_3d[i][0]
    ell_2=coordinates_tuple_3d[i][2]-coordinates_tuple_3d[i][0]
    ell_1_mag=compute_vector_magnitude(ell_1)
    ell_2_mag=compute_vector_magnitude(ell_2)
    area_vector=np.cross(ell_1,ell_2)
    area_vector_mag=compute_vector_magnitude(area_vector)

    # area vector and rotation_vector should point in same direction 

    if m.isclose(np.dot(compute_unit_vector(area_vector),compute_unit_vector(Rotation_vector[:,i])),1):
        
        x=area_vector[0]
        y=area_vector[1]
        z=area_vector[2]
        
        # radial coord
        spherical_coordinates_area_vector[i,0]=np.sqrt((x**2)+(y**2)+(z**2))
        # theta coord azimuth
        spherical_coordinates_area_vector[i,1]=np.sign(y)*np.arccos(x/(np.sqrt((x**2)+(y**2)))) 
        # phi coord incline 
        spherical_coordinates_area_vector[i,2]=np.arccos((z)/np.sqrt((x**2)+(y**2)+(z**2)))
    else:
        print("area vector and rotation vector do not align")
        print("index=",i)
        print("dot product=",np.dot(compute_unit_vector(area_vector),compute_unit_vector(Rotation_vector[:,i])))
        break
        


    # x=area_vector[0]
    # y=area_vector[1]
    # z=area_vector[2]
    
    # # radial coord
    # spherical_coordinates_area_vector[i,0]=np.sqrt((x**2)+(y**2)+(z**2))
    # # theta coord azimuth
    # spherical_coordinates_area_vector[i,1]=np.sign(y)*np.arccos(x/(np.sqrt((x**2)+(y**2)))) 
    # # phi coord incline 
    # spherical_coordinates_area_vector[i,2]=np.arccos((z)/np.sqrt((x**2)+(y**2)+(z**2)))


    
#%% plot theta histogram 

pi_theta_ticks=[ -np.pi, -np.pi/2, 0, np.pi/2,np.pi]
pi_theta_tick_labels=['-π','-π/2','0', 'π/2', 'π'] 
plt.hist((spherical_coordinates_area_vector[:,1]),density=True, bins=bin_count)
plt.xticks(pi_theta_ticks, pi_theta_tick_labels)
plt.xlabel('$\\theta$')
plt.show()

pi_phi_ticks=[ 0,np.pi/4, np.pi/2,3*np.pi/4,np.pi]
pi_phi_tick_labels=[ '0','π/4', 'π/2','3π/4' ,'π']
plt.hist(spherical_coordinates_area_vector[:,2],density=True, bins=bin_count)
plt.xticks(pi_phi_ticks,pi_phi_tick_labels)
plt.xlabel('$\phi$')
plt.show()






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
#abs_path_2_lammps_exec='/Users/luke_dev/Documents/lammps_hirotori/build_serial_maxbook/lmp'
# abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5'
#abs_path_2_lammps_script='/Users/luke_dev/Documents/LSC/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5'

# for running on myriad 
abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps_hirotori/build_serial/lmp'
abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5'


Path_2_generic='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'



extra_code='module unload mpi compilers gcc-libs \n module load beta-modules \n module load gcc-libs/10.2.0 \n module load compilers/intel/2022.2 \n module load mpi/intel/2019/update6/intel \n module load hdf/5-1.10.6/gnu-10.2.0'
wd_path='/home/ucahlrl/Scratch/output/'
num_task_req=''
data_transfer_instructions=''
SRD_MD_ratio_ = 10
VP_ave_freq=10000
md_timestep=0.005071624521210362
collision_time_negative_bar=0.05071624521210362
#erate= np.array([0.0008,0.001,0.002,0.005,0.01,0.1])
erate= np.array([0.03,0.0275,0.025,0.0225,0.02,0.0175,0.015,0.0125,0.01,0.0075])#0.0075,0.005,0.0025,0.001,0.0001]) 
#erate=np.array([0.03,0.0275,0.025])
#erate= np.array([0.0075,0.005,0.0025]) longer runs which need checkpointing
dump_freq_=10000
thermo_freq=10000#dump_freq_
i_=0
j_=number_of_points
fluid_name='langevinrun'
bending_stiffness=np.array([10000]) # original 50,100,200,400
#internal_stiffness=np.array([60,80,100]) # 20 does nothing 
internal_stiffness=np.array([60])
#no_timestep_=np.array([2000000]) 
# just for dist test 
#no_timestep_=np.array([1000]) 
#phantom_mass=np.array([0.01])
equilibrium_triangle_side_length=3
tempdir_req='1G'
ram_requirement='2G'
wall_time='00:30:00'
damp=0.05

np_req=np.array([1,1,1,1,1,1,1,1,1,1]).astype('str')
initial_temp=0.81
inital_temp_multipl=np.array([1, 1, 1, 1, 1,
       1, 1, 1,1, 1])

input_temp=initial_temp*inital_temp_multipl

realisation_index_=[1,2,3]
realisation_index_=np.arange(0,1000,1)
timestep_multiplier=0.05

def compute_timesteps_for_strain(total_strain,erate,md_timestep,timestep_multiplier):
      no_timestep_=(np.round((total_strain/(erate*md_timestep*timestep_multiplier)),decimals=-3)).astype('int')

      return no_timestep_

# for MYRIAD run
total_strain=600
no_timestep_=compute_timesteps_for_strain(total_strain,erate,md_timestep,timestep_multiplier)




# # for bug test
# number_of_restarts_per_run=np.array([1,1,1,1,1])
# no_timestep_=np.array([2000,2000,2000,2000,2000])
# np_req=np.array([14,14,14,14,14]).astype('str')






# for 8,16,32,36
var_choice_1=erate
var_choice_2=internal_stiffness
# individual shear rates 
def sim_file_prod_flat_elastic_initial_MYRIAD(input_temp,coordinates_tuple_3d,erate,equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_eqts_'+str(equilibrium_triangle_side_length)+'_realisations_'+str(j_)+'_box_size_'+str(box_size_bar[box_size_index])+'_damp_'+str(damp)+'_bendstiff_'+str(bending_stiffness[0])+'_intstiff_'+str(internal_stiffness[0])+'_'+str(internal_stiffness[-1])+'_erate_'+str(erate)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    run_code_list=[]
    # test to check consistency of cores request 

    

    #for n in range(0,np_req.size):
        #or now just use one realisation 
    #for k in range(0,var_choice_1.size):    
    #for k in range(0,1):   
    for m in range(0,var_choice_2.size): 
                #for m in range(0,1):  
                        for j in range(i_,j_):
                            param_set_code=str(np.random.randint(0, 1000000))
                            simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_realisation_'+str(j)+'_Bk_'+str(bending_stiffness[0])+'_np_'+str(np_req[k])+'_no_timesteps_'+str(no_timestep_[k])+'_intstiff_'+str(var_choice_2[m])+'_eqsl_'+str(equilibrium_triangle_side_length)+'_erate_'+str(erate)+'_'
                            run_code=''
                           
                            #print(no_SRD)
                            box_size = str(box_size_bar[box_size_index])
                            timestep_input= str(md_timestep)
                            # number of chunks to use for VP averaging
                            SRD_MD_ratio=str(int(SRD_MD_ratio_))
                            lamda= str(collision_time_negative_bar)
                            dump_freq=str(10000)
                            thermo_freq = str(10000)
                            no_timesteps = str(no_timestep_[n])
                            rand_int =str(np.random.randint(0, 1000000))
                            rand_int_1 =str( np.random.randint(0, 1000000))
                            rand_int_2 =str(np.random.randint(0, 1000000))
                            rand_int_3=str(np.random.randint(0,1000000))
                            num_proc=str(np_req[k])
                          

                            stokes_bead_1=coordinates_tuple_3d[j][0]
                            stokes_bead_2=coordinates_tuple_3d[j][1]
                            stokes_bead_3=coordinates_tuple_3d[j][2]
                            phantom_bead_1=coordinates_tuple_3d[j][3] 
                            phantom_bead_2=coordinates_tuple_3d[j][4]
                            phantom_bead_3=coordinates_tuple_3d[j][5]

                                        
                            run_code_individual =abs_path_2_lammps_exec+' -var temp '+str(input_temp)+' -var damp '+str(damp)+' -var erate_in '+str(erate)+' -var equilirbium_triangle_side_length '+str(equilibrium_triangle_side_length)+\
                            ' -var bead_1_x_position '+str(stokes_bead_1[0])+' -var bead_1_y_position '+str(stokes_bead_1[1])+' -var bead_1_z_position '+str(stokes_bead_1[2])+' -var bead_2_x_position '+str(stokes_bead_2[0])+\
                            ' -var bead_2_y_position '+str(stokes_bead_2[1])+' -var bead_2_z_position '+str(stokes_bead_2[2])+' -var bead_3_x_position '+str(stokes_bead_3[0])+' -var bead_3_y_position '+str(stokes_bead_3[1])+\
                            ' -var bead_3_z_position '+str(stokes_bead_3[2])+' -var bead_1p_x_position '+str(phantom_bead_1[0])+' -var bead_1p_y_position '+str(phantom_bead_1[1])+' -var bead_1p_z_position '+str(phantom_bead_1[2])+\
                            ' -var bead_2p_x_position '+str(phantom_bead_2[0])+' -var bead_2p_y_position '+str(phantom_bead_2[1])+' -var bead_2p_z_position '+str(phantom_bead_2[2])+' -var bead_3p_x_position '+str(phantom_bead_3[0])+\
                            ' -var bead_3p_y_position '+str(phantom_bead_3[1])+' -var bead_3p_z_position '+str(phantom_bead_3[2])+\
                            ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[m])+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+\
                            ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[j])+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
                            ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
                            ' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -in '+abs_path_2_lammps_script+' \n '  #>> '+prod_run_file_name+' & \n'
                            
                            run_code_list.append(run_code_individual)
                            run_code=run_code +run_code_individual

                            
                            run_code = run_code[:-2]

                            py2bash_launch_overwriter.py2bash_launch_overwriter_serial(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,num_proc,wd_path,extra_code,run_code,data_transfer_instructions)
    return run_code_list, sim_batchcode

def sim_file_prod_flat_elastic_MYRIAD_all_erate_one_file(input_temp,coordinates_tuple_3d,erate,equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_eqts_'+str(equilibrium_triangle_side_length)+'_realisations_'+str(j_)+'_box_size_'+str(box_size_bar[box_size_index])+'_damp_'+str(damp)+'_bendstiff_'+str(bending_stiffness[0])+'_intstiff_'+str(internal_stiffness[0])+'_'+str(internal_stiffness[-1])+'_erate_'+str(erate[0])+'_'+str(erate[-1])+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    run_code_list=[]
    # test to check consistency of cores request 

    

    #for n in range(0,np_req.size):
        #or now just use one realisation 
    for k in range(0,var_choice_1.size):    
    #for k in range(0,1):   
        for m in range(0,var_choice_2.size): 
                    #for m in range(0,1):  
                            for j in range(i_,j_):
                                param_set_code=str(np.random.randint(0, 1000000))
                                simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_realisation_'+str(j)+'_Bk_'+str(bending_stiffness[0])+'_np_'+str(np_req[k])+'_no_timesteps_'+str(no_timestep_[k])+'_intstiff_'+str(var_choice_2[m])+'_eqsl_'+str(equilibrium_triangle_side_length)+'_erate_'+str(erate[k])+'_'
                                run_code=''
                            
                                #print(no_SRD)
                                box_size = str(box_size_bar[box_size_index])
                                timestep_input= str(md_timestep)
                                # number of chunks to use for VP averaging
                                SRD_MD_ratio=str(int(SRD_MD_ratio_))
                                lamda= str(collision_time_negative_bar)
                                dump_freq=str(10000)
                                thermo_freq = str(10000)
                                no_timesteps = str(no_timestep_[k])
                                rand_int =str(np.random.randint(0, 1000000))
                                rand_int_1 =str( np.random.randint(0, 1000000))
                                rand_int_2 =str(np.random.randint(0, 1000000))
                                rand_int_3=str(np.random.randint(0,1000000))
                                num_proc=str(np_req[k])
                            

                                stokes_bead_1=coordinates_tuple_3d[j][0]
                                stokes_bead_2=coordinates_tuple_3d[j][1]
                                stokes_bead_3=coordinates_tuple_3d[j][2]
                                phantom_bead_1=coordinates_tuple_3d[j][3] 
                                phantom_bead_2=coordinates_tuple_3d[j][4]
                                phantom_bead_3=coordinates_tuple_3d[j][5]

                                            
                                run_code_individual =abs_path_2_lammps_exec+' -var temp '+str(input_temp[k])+' -var damp '+str(damp)+' -var erate_in '+str(erate[k])+' -var equilirbium_triangle_side_length '+str(equilibrium_triangle_side_length)+\
                                ' -var bead_1_x_position '+str(stokes_bead_1[0])+' -var bead_1_y_position '+str(stokes_bead_1[1])+' -var bead_1_z_position '+str(stokes_bead_1[2])+' -var bead_2_x_position '+str(stokes_bead_2[0])+\
                                ' -var bead_2_y_position '+str(stokes_bead_2[1])+' -var bead_2_z_position '+str(stokes_bead_2[2])+' -var bead_3_x_position '+str(stokes_bead_3[0])+' -var bead_3_y_position '+str(stokes_bead_3[1])+\
                                ' -var bead_3_z_position '+str(stokes_bead_3[2])+' -var bead_1p_x_position '+str(phantom_bead_1[0])+' -var bead_1p_y_position '+str(phantom_bead_1[1])+' -var bead_1p_z_position '+str(phantom_bead_1[2])+\
                                ' -var bead_2p_x_position '+str(phantom_bead_2[0])+' -var bead_2p_y_position '+str(phantom_bead_2[1])+' -var bead_2p_z_position '+str(phantom_bead_2[2])+' -var bead_3p_x_position '+str(phantom_bead_3[0])+\
                                ' -var bead_3p_y_position '+str(phantom_bead_3[1])+' -var bead_3p_z_position '+str(phantom_bead_3[2])+\
                                ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[m])+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+\
                                ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[j])+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
                                ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
                                ' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -in '+abs_path_2_lammps_script+' \n '  #>> '+prod_run_file_name+' & \n'
                                
                                run_code_list.append(run_code_individual)
                                run_code=run_code +run_code_individual

                                
                                run_code = run_code[:-2]

                                py2bash_launch_overwriter.py2bash_launch_overwriter_serial(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,num_proc,wd_path,extra_code,run_code,data_transfer_instructions)
    return run_code_list, sim_batchcode


# %%
os.chdir("/Users/luke_dev/Documents/simulation_test_folder/")
k=0
#os.chdir("/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder/bug_test_spring_force/")
for n in range(erate.size): 
    os.chdir("/Users/luke_dev/Documents/simulation_test_folder/")
    MyFile=open('K_'+str(internal_stiffness[0])+'_damp_'+str(damp)+'_erate_'+str(erate[n])+'.sh','w')
    outputs= sim_file_prod_flat_elastic_initial_MYRIAD(input_temp[n],
                                                    coordinates_tuple_3d,
                                                    erate[n],
                                                    equilibrium_triangle_side_length,
                                                    var_choice_1,
                                                    var_choice_2,
                                                    internal_stiffness,
                                                    data_transfer_instructions,
                                                    extra_code,
                                                    wd_path,
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
                                                    md_timestep,
                                                    i_,j_,box_size_bar,
                                                    box_size_index,
                                                    Path_2_shell_scirpts,
                                                    Path_2_generic,
                                                    fluid_name)
    run_code_list=outputs[0]
    sim_batchcode=outputs[1]
    for element in run_code_list:
        MyFile.write(element)
        MyFile.write('\n')
    MyFile.close()


# %%
# all in one file for myriad

sim_file_prod_flat_elastic_MYRIAD_all_erate_one_file(input_temp,
                                                     coordinates_tuple_3d,
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
                                                     md_timestep,
                                                     i_,
                                                     j_,
                                                     box_size_bar,
                                                     box_size_index,
                                                     Path_2_shell_scirpts,
                                                     Path_2_generic,
                                                     fluid_name)
    

# %%
