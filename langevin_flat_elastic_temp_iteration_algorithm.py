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


number_of_points=10 # for bug test , normally 30
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


Path_2_shell_scirpts='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'

# for running on my computer 
abs_path_2_lammps_exec='/Users/luke_dev/Documents/lammps_hirotori/build_serial_maxbook/lmp'
# abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5'
abs_path_2_lammps_script='/Users/luke_dev/Documents/LSC/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5'
abs_path_2_lammps_script='/Users/luke_dev/Documents/LSC/in.langevin_with_hookean_flat_elastic_particle_only_dump_hdf5_eq_b4_shear'


Path_2_generic='/Users/luke_dev/Documents/Shell_scripts_for_MYRIAD'


num_task_req=''
data_transfer_instructions=''
SRD_MD_ratio_ = 10
VP_ave_freq=10000
md_timestep=0.005071624521210362
collision_time_negative_bar=0.05071624521210362


#erate= np.array([0.0075,0.005,0.0025]) longer runs which need checkpointing
erate=np.array([1,0.9,0.7,0.5,0.2,0.1,0.09,0.08,
                0.07,0.06,0.05,0.04,
                0.03,0.0275,0.025,0.0225,
                0.02,0.0175,0.015,0.0125,
                0.01,0.0075,0.005,0.0025,
                0.001,0.00075,0.0005])
# erate=np.array([100,50,25,10,5])
dump_freq=str(100)
thermo_freq=str(100)#dump_freq_
i_=0
j_=number_of_points
fluid_name='langevinrun'
bending_stiffness=np.array([10000]) # original 50,100,200,400


damp_init=0.03633 # based on H20 as reference fluid 


internal_stiffness=np.array([1500,2000,2500])




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
total_strain=50
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
#[0.18999999999999928,
#  0.4299999999999995,
#  0.5299999999999996,
#  0.5799999999999996,
#  0.6499999999999997,
# 0.6599999999999997,
#  0.6599999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
# 0.6699999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
# 0.6699999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
#  0.6699999999999997,
# 0.6699999999999997,
# 0.6699999999999997, 
# 0.6699999999999997
# 0.6699999999999997,
# 0.6699999999999997
#0.6699999999999997
#0.6699999999999997],

count=0
tolerance=0.01
increment=0.01
input_temp=1
desired_temp=1
from log2numpy import *
thermo_vars='         KinEng         PotEng        c_myTemp        c_bias         TotEng    '
filepath='/Users/luke_dev/Documents/simulation_test_folder/temp_iteration_test'
os.chdir(filepath)
temp_comparison=0.5
new_temp=[]
for n in range(5,10):
    input_temp=1
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
        box_size = str(box_size_bar[box_size_index])
        timestep_input= str(md_timestep)
        # number of chunks to use for VP averaging
        SRD_MD_ratio=str(int(SRD_MD_ratio_))
        lamda= str(collision_time_negative_bar)

        no_timesteps = str(no_timestep_[n])
        rand_int =str(np.random.randint(0, 1000000))
        rand_int_1 =str( np.random.randint(0, 1000000))
        rand_int_2 =str(np.random.randint(0, 1000000))
        rand_int_3=str(np.random.randint(0,1000000))



        stokes_bead_1=coordinates_tuple_3d[j][0]
        stokes_bead_2=coordinates_tuple_3d[j][1]
        stokes_bead_3=coordinates_tuple_3d[j][2]
        phantom_bead_1=coordinates_tuple_3d[j][3] 
        phantom_bead_2=coordinates_tuple_3d[j][4]
        phantom_bead_3=coordinates_tuple_3d[j][5]

                    
        run_code_individual =abs_path_2_lammps_exec+' -var temp '+str(input_temp)+' -var damp '\
            +str(damp_init)+' -var erate_in '+str(erate_in)+' -var equilirbium_triangle_side_length '\
                +str(equilibrium_triangle_side_length)+\
        ' -var bead_1_x_position '+str(stokes_bead_1[0])+' -var bead_1_y_position '+str(stokes_bead_1[1])+\
            ' -var bead_1_z_position '+str(stokes_bead_1[2])+' -var bead_2_x_position '+str(stokes_bead_2[0])+\
        ' -var bead_2_y_position '+str(stokes_bead_2[1])+' -var bead_2_z_position '+str(stokes_bead_2[2])+\
            ' -var bead_3_x_position '+str(stokes_bead_3[0])+' -var bead_3_y_position '+str(stokes_bead_3[1])+\
        ' -var bead_3_z_position '+str(stokes_bead_3[2])+' -var bead_1p_x_position '+str(phantom_bead_1[0])+\
            ' -var bead_1p_y_position '+str(phantom_bead_1[1])+' -var bead_1p_z_position '+str(phantom_bead_1[2])+\
        ' -var bead_2p_x_position '+str(phantom_bead_2[0])+' -var bead_2p_y_position '+str(phantom_bead_2[1])+\
            ' -var bead_2p_z_position '+str(phantom_bead_2[2])+' -var bead_3p_x_position '+str(phantom_bead_3[0])+\
        ' -var bead_3p_y_position '+str(phantom_bead_3[1])+' -var bead_3p_z_position '+str(phantom_bead_3[2])+\
        ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[l])+\
            ' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(999)+\
        ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[j])+\
            ' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
        ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var box_size '+box_size+\
            ' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
        ' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+\
            no_timesteps+' -in '+abs_path_2_lammps_script+' \n '  #>> '+prod_run_file_name+' & \n'

        os.system(run_code_individual)

        logfile_name="log.langevinrun_no"+str(999)+"_hookean_flat_elastic_"+rand_int+\
            "_0_23_0.03633_0.005071624521210362_100_100_"+no_timesteps+"_0.2_gdot_"+str(erate_in)+\
            "_BK_10000_K_"+str(internal_stiffness[l])
        log_file_data=log2numpy_reader(logfile_name,
                                    filepath,
                                    thermo_vars)

        temp_out=np.mean(log_file_data[5:,4])
        temp_comparison=desired_temp-temp_out
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
            new_temp.append(input_temp)


     
     
         

     
#log.${fluid_name}_no${sim_batchcode}_hookean_flat_elastic_${rand_int}_${realisation_index}_${box_size}_${damp}_${timestep_input}_${dump_freq}_${thermo_freq}_${no_timesteps}_${timestep_multiplier}_gdot_${erate_in}_BK_${angle_stiff}_K_${spring_stiffness}




    



# %% prodicing files in simulation test folder

os.chdir("/Users/luke_dev/Documents/simulation_test_folder/")
kdamp=int(internal_stiffness_init*damp_init)
path="/Users/luke_dev/Documents/simulation_test_folder/low_damp_test"
# os.mkdir(path)
os.chdir(path)

k=0
#os.chdir("/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder/bug_test_spring_force/")



# %%
# all in one file for myriad

        

# %%
