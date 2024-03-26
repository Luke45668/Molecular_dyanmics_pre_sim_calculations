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
print(collision_time_negative_bar)
print("collision ratio",collision_ratio)
print("box size",box_size_bar)
print("srd count",srd_count)
print("estimated run time in hours", esimtated_run_time_mins/60)
print("estimated equilibration steps", estimated_equilibration_steps)

#%% Producing random dist of orientations 
# produce N random numbers between 0,1 

#NOTE using mathematics convention for the angles
# r 
# theta is the equator 
#phi is the incline


number_of_points=30
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

box_size_index=7

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


#%% KATHLEEN paths 
##################


Path_2_shell_scirpts='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_KATHLEEN'
abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps-23Jun2022_with_SRD_pol/build_kathleen/lmp_mpi'
abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/pure_MPCD_with_hookean_flat_elastic.file'
Path_2_generic='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_KATHLEEN'
data_transfer_instructions=''
extra_code=''
wd_path='/home/ucahlrl/Scratch/output/'
num_task_req=''
tempdir_req=''
##################
### Variables settings
# calculating possible SRD_MD ratios and MD timesteps 
SRD_MD_ratio_ = np.array([1000])
md_timestep=collision_time_negative_bar/SRD_MD_ratio_
no_timestep_=(10000000*(SRD_MD_ratio_/SRD_MD_ratio_[0])).astype('int')
print(md_timestep)
thermo_freq=1000
VP_ave_freq=10000
i_=0
j_=3
Path_2_generic='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_KATHLEEN'
hypthread='2'
fluid_name='flatelastictest'
max_cores= 64
swap_number=1
swap_rate=np.array([15])
bending_stiffness=np.array([100,200,300]) # original 50,100,200,400
internal_stiffness=np.array([10])
phantom_mass=np.array([0.01])
equilibrium_triangle_side_length=3
box_size_index=4
realisation_index_=[1,2,3]

ram_requirement='2G'
wall_time='24:00:00'
np_req=str(max_cores)
phi_= str(phi[box_size_index])
 # for one KATHLEEN node 
if (int(np_req)) > max_cores:
      print("Too many cores requested")
      breakpoint()
else:
      print("Core request satisfactory, producing simulation submission script ")

var_choice_1=bending_stiffness
var_choice_2=SRD_MD_ratio_

def sim_file_prod_flat_elastic_individual_kathleen(equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,hypthread,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,swap_number,i_,j_,swap_rate,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_eqts_'+str(equilibrium_triangle_side_length)+'_fluid_visc_'+str(nu_bar)+'_box_size_'+str(box_size_bar[box_size_index])+'_swap_rate_'+str(swap_rate[0])+'_intstiff_'+str(internal_stiffness[0])+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    # test to check consistency of cores request 
    
    
    for j in range(i_,j_): #or now just use one realisation 
        for k in range(0,var_choice_1.size):       
                for m in range(0,var_choice_2.size):#range(0,1):  
                        param_set_code=str(np.random.randint(0, 1000000))
                        simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_realisation_'+str(j)+'_Bk_'+str(var_choice_1[k])+'_timestep_'+str(md_timestep[m])+'_no_timesteps_'+str(no_timestep_[m])+'_SRDMDratio_'+str(var_choice_2[m])+'_eqts_'+str(equilibrium_triangle_side_length)+'_'
                        run_code=''
                        no_SRD=str(int(srd_count[box_size_index])) 
                        #print(no_SRD)
                        box_size = str(box_size_bar[box_size_index])
                        timestep_input= str(md_timestep[m])
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=str(int(SRD_MD_ratio_[m]))
                        lamda= str(collision_time_negative_bar)
                        dump_freq=str(int(SRD_MD_ratio_[m]/2)) # 
                        thermo_freq = str(thermo_freq) # 
                        no_timesteps = str(no_timestep_[m])
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))
                        num_proc=str(np_req)
                        


                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var equilirbium_triangle_side_length '+str(equilibrium_triangle_side_length)+' -var angle_stiff '+str(bending_stiffness[k])+' -var spring_stiffness '+str(internal_stiffness[0])+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[0])+' -var swap_number '+str(swap_number)+' -var VP_ave_freq '+str(VP_ave_freq)+' -var chunk '+str(chunk)+' -var realisation_index '+str(realisation_index_[j])+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -in '+abs_path_2_lammps_script+' \n &'  #>> '+prod_run_file_name+' & \n'
                        #print(run_code_individual)
                        run_code=run_code +run_code_individual


                        run_code = run_code[:-3]

                        py2bash_launch_overwriter_kathleen.py2bash_launch_overwriter(hypthread,Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)


# %%
sim_file_prod_flat_elastic_individual_kathleen(equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,hypthread,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,swap_number,i_,j_,swap_rate,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name)


# %% producing shell scripts for MYRIAD with rotations 
# on my computer 
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_imac_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_macbook_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_macbook_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_imac_h5md/lmp_mpi'
#abs_path_2_lammps_script='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/LSC/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5'
#abs_path_2_lammps_script='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/LSC/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5_chain'


Path_2_shell_scirpts='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD'

abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps_hirotori/build_MYRIAD/lmp_mpi'
# abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5'
abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5_chain'
Path_2_generic='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD'



extra_code='module unload mpi compilers gcc-libs \n module load beta-modules \n module load gcc-libs/10.2.0 \n module load compilers/intel/2022.2 \n module load mpi/intel/2019/update6/intel \n module load hdf/5-1.12.3-impi/intel-2022'
wd_path='/home/ucahlrl/Scratch/output/'
num_task_req=''
data_transfer_instructions=''
SRD_MD_ratio_ = 10
VP_ave_freq=10000
md_timestep=collision_time_negative_bar/SRD_MD_ratio_
#erate= np.array([0.0008,0.001,0.002,0.005,0.01,0.1])
erate= np.array([0.02,0.0175,0.015,0.0125,0.01]) 
#erate= np.array([0.0075,0.005,0.0025]) longer runs which need checkpointing
dump_freq_=10
thermo_freq=1000
i_=0
j_=number_of_points
fluid_name='feshortruns'
bending_stiffness=np.array([10000]) # original 50,100,200,400
internal_stiffness=np.array([60,80,100]) # 20 does nothing 
#no_timestep_=np.array([2000000]) 
# just for dist test 
#no_timestep_=np.array([1000]) 
#phantom_mass=np.array([0.01])
equilibrium_triangle_side_length=3
tempdir_req='1G'
ram_requirement='8G'
wall_time='48:00:00'
np_req=np.array([16,16,16,16,16]).astype('str')

phi_= str(phi[box_size_index])
realisation_index_=[1,2,3]
realisation_index_=np.arange(0,1000,1)
timestep_multiplier=0.05

def compute_timesteps_for_strain(total_strain,erate,md_timestep,timestep_multiplier):
      no_timestep_=(np.round((total_strain/(erate*md_timestep*timestep_multiplier)),decimals=-3)).astype('int')

      return no_timestep_

# for MYRIAD run
total_strain=30
no_timestep_=compute_timesteps_for_strain(total_strain,erate,md_timestep,timestep_multiplier)
number_of_restarts_per_run=np.array([8,8,8,8,8])
restart_frequency=(no_timestep_/number_of_restarts_per_run).astype('int')
end_run_strain=()
for i in range(erate.size):
    end_run_strain=end_run_strain+ (np.arange(1,number_of_restarts_per_run[i]+1)*(total_strain/number_of_restarts_per_run[i]), )
# for test run 
#no_timestep_=np.array([10000,10000,10000,10000,10000,10000])
# to guide processor choice
max_number_of_steps_possible=np.floor(np.array([ 8627174.4, 18864726.4 , 43018188.8 ,47761926.40000001]))
# for 8,16,32,36
var_choice_1=erate
var_choice_2=internal_stiffness

def sim_file_prod_flat_elastic_initial_MYRIAD(restart_frequency,coordinates_tuple_3d,erate,equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_eqts_'+str(equilibrium_triangle_side_length)+'_realisations_'+str(j_)+'_box_size_'+str(box_size_bar[box_size_index])+'_bendstiff_'+str(bending_stiffness[0])+'_intstiff_'+str(internal_stiffness[0])+'_'+str(internal_stiffness[-1])+'_'+ META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    run_code_list=[]
    # test to check consistency of cores request 

    

    #for n in range(0,np_req.size):
        #or now just use one realisation 
    for k in range(0,var_choice_1.size):    
        #for k in range(0,3):   
                for m in range(0,var_choice_2.size): 
                #for m in range(0,1):  
                        for j in range(i_,j_):
                            param_set_code=str(np.random.randint(0, 1000000))
                            simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_realisation_'+str(j)+'_Bk_'+str(bending_stiffness[0])+'_np_'+str(np_req[k])+'_no_timesteps_'+str(no_timestep_[k])+'_intstiff_'+str(var_choice_2[m])+'_eqsl_'+str(equilibrium_triangle_side_length)+'_erate_'+str(erate[k])+'_'
                            run_code=''
                            no_SRD=str(int(srd_count[box_size_index])) 
                            #print(no_SRD)
                            box_size = str(box_size_bar[box_size_index])
                            timestep_input= str(md_timestep)
                            # number of chunks to use for VP averaging
                            SRD_MD_ratio=str(int(SRD_MD_ratio_))
                            lamda= str(collision_time_negative_bar)
                            dump_freq=str(int(SRD_MD_ratio_)) # 
                            thermo_freq = str(thermo_freq) # 
                            no_timesteps = str(no_timestep_[k])
                            rand_int =str(np.random.randint(0, 1000000))
                            rand_int_1 =str( np.random.randint(0, 1000000))
                            rand_int_2 =str(np.random.randint(0, 1000000))
                            rand_int_3=str(np.random.randint(0,1000000))
                            num_proc=str(np_req[k])
                            restart_freq=restart_frequency[k]
                            end_strain=total_strain

                            stokes_bead_1=coordinates_tuple_3d[j][0]
                            stokes_bead_2=coordinates_tuple_3d[j][1]
                            stokes_bead_3=coordinates_tuple_3d[j][2]
                            phantom_bead_1=coordinates_tuple_3d[j][3] 
                            phantom_bead_2=coordinates_tuple_3d[j][4]
                            phantom_bead_3=coordinates_tuple_3d[j][5]

                                        
                            run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var end_run_strain '+str(end_strain)+' -var restart_freq '+str(restart_freq)+' -var erate_in '+str(erate[k])+' -var equilirbium_triangle_side_length '+str(equilibrium_triangle_side_length)+\
                            ' -var bead_1_x_position '+str(stokes_bead_1[0])+' -var bead_1_y_position '+str(stokes_bead_1[1])+' -var bead_1_z_position '+str(stokes_bead_1[2])+' -var bead_2_x_position '+str(stokes_bead_2[0])+\
                            ' -var bead_2_y_position '+str(stokes_bead_2[1])+' -var bead_2_z_position '+str(stokes_bead_2[2])+' -var bead_3_x_position '+str(stokes_bead_3[0])+' -var bead_3_y_position '+str(stokes_bead_3[1])+\
                            ' -var bead_3_z_position '+str(stokes_bead_3[2])+' -var bead_1p_x_position '+str(phantom_bead_1[0])+' -var bead_1p_y_position '+str(phantom_bead_1[1])+' -var bead_1p_z_position '+str(phantom_bead_1[2])+\
                            ' -var bead_2p_x_position '+str(phantom_bead_2[0])+' -var bead_2p_y_position '+str(phantom_bead_2[1])+' -var bead_2p_z_position '+str(phantom_bead_2[2])+' -var bead_3p_x_position '+str(phantom_bead_3[0])+\
                            ' -var bead_3p_y_position '+str(phantom_bead_3[1])+' -var bead_3p_z_position '+str(phantom_bead_3[2])+\
                            ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[m])+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+\
                            ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[j])+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
                            ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var no_SRD '+no_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
                            ' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -in '+abs_path_2_lammps_script+' \n &'  #>> '+prod_run_file_name+' & \n'
                            
                            run_code_list.append(run_code_individual)
                            run_code=run_code +run_code_individual

                            
                            run_code = run_code[:-3]

                            py2bash_launch_overwriter.py2bash_launch_overwriter(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,num_proc,wd_path,extra_code,run_code,data_transfer_instructions)
    return run_code_list, sim_batchcode



# %%
outputs= sim_file_prod_flat_elastic_initial_MYRIAD(restart_frequency,coordinates_tuple_3d,erate,equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name)
run_code_list=outputs[0]
sim_batchcode=outputs[1]

# %%
#os.mkdir("/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder/test_run_flat_elastic_100_realisations/")
os.chdir("/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/")
MyFile=open('initial_chain_test.sh','w')

for element in run_code_list:
     MyFile.write(element)
     MyFile.write('\n')
MyFile.close()

# %% producing subsequent restart files
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_imac_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_macbook_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_macbook_h5md/lmp_mpi'
abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_imac_h5md/lmp_mpi'


abs_path_2_lammps_script='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/LSC/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5_restart_only'



#Path_2_shell_scirpts='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD'

#abs_path_2_lammps_exec='/home/ucahlrl/simulation_run_folder/lammps_hirotori/build_MYRIAD/lmp_mpi'
#abs_path_2_lammps_script='/home/ucahlrl/simulation_run_folder/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5_restart_only'
Path_2_generic='/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD'




strain_start_point=30
restart_count=1 
previous_end_run_strain=()
for i in range(erate.size):
    previous_end_run_strain=previous_end_run_strain+ (strain_start_point+np.arange(1,number_of_restarts_per_run[i]+1)*(total_strain/number_of_restarts_per_run[i]), )




def sim_file_prod_flat_elastic_restarts_MYRIAD(restart_count,sim_batchcode,strain_start_point,restart_frequency,coordinates_tuple_3d,erate,equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_strain_start_'+str(strain_start_point)+'_realisations_'+str(j_)+'_box_size_'+str(box_size_bar[box_size_index])+'_bendstiff_'+str(bending_stiffness[0])+'_intstiff_'+str(internal_stiffness[0])+'_'+str(internal_stiffness[-1])+'_'+ META_DATA
    os.mkdir(simulation_batch_folder)
    #sim_batchcode=str(np.random.randint(0, 1000000))
    run_code_list=[]
    # test to check consistency of cores request 

    

    #for n in range(0,np_req.size):
        #or now just use one realisation 
    for k in range(0,var_choice_1.size):    
        #for k in range(0,3):   
                for m in range(0,var_choice_2.size): 
                #for m in range(0,1):  
                        for j in range(i_,j_):
                            param_set_code=str(np.random.randint(0, 1000000))
                            simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_realisation_'+str(j)+'_Bk_'+str(bending_stiffness[0])+'_np_'+str(np_req[k])+'_no_timesteps_'+str(no_timestep_[k])+'_intstiff_'+str(var_choice_2[m])+'_eqsl_'+str(equilibrium_triangle_side_length)+'_erate_'+str(erate[k])+'_'
                            run_code=''
                            no_SRD=str(int(srd_count[box_size_index])) 
                            #print(no_SRD)
                            box_size = str(box_size_bar[box_size_index])
                            timestep_input= str(md_timestep)
                            # number of chunks to use for VP averaging
                            SRD_MD_ratio=str(int(SRD_MD_ratio_))
                            lamda= str(collision_time_negative_bar)
                            dump_freq=str(int(SRD_MD_ratio_)) # 
                            thermo_freq = str(thermo_freq) # 
                            no_timesteps = str(no_timestep_[k])
                            rand_int =str(np.random.randint(0, 1000000))
                            rand_int_1 =str( np.random.randint(0, 1000000))
                            rand_int_2 =str(np.random.randint(0, 1000000))
                            rand_int_3=str(np.random.randint(0,1000000))
                            num_proc=str(np_req[k])
                            restart_freq=restart_frequency[k]
                            previous_end_strain=strain_start_point
                            end_strain=strain_start_point+total_strain

                            stokes_bead_1=coordinates_tuple_3d[j][0]
                            stokes_bead_2=coordinates_tuple_3d[j][1]
                            stokes_bead_3=coordinates_tuple_3d[j][2]
                            phantom_bead_1=coordinates_tuple_3d[j][3] 
                            phantom_bead_2=coordinates_tuple_3d[j][4]
                            phantom_bead_3=coordinates_tuple_3d[j][5]

                                        
                            run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var restart_count '+str(restart_count)+' -var previous_end_strain '+str(strain_start_point)+' -var end_run_strain '+str(end_strain)+\
                            '  -var restart_freq '+str(restart_freq)+' -var erate_in '+str(erate[k])+' -var equilirbium_triangle_side_length '+str(equilibrium_triangle_side_length)+\
                            ' -var bead_1_x_position '+str(stokes_bead_1[0])+' -var bead_1_y_position '+str(stokes_bead_1[1])+' -var bead_1_z_position '+str(stokes_bead_1[2])+' -var bead_2_x_position '+str(stokes_bead_2[0])+\
                            ' -var bead_2_y_position '+str(stokes_bead_2[1])+' -var bead_2_z_position '+str(stokes_bead_2[2])+' -var bead_3_x_position '+str(stokes_bead_3[0])+' -var bead_3_y_position '+str(stokes_bead_3[1])+\
                            ' -var bead_3_z_position '+str(stokes_bead_3[2])+' -var bead_1p_x_position '+str(phantom_bead_1[0])+' -var bead_1p_y_position '+str(phantom_bead_1[1])+' -var bead_1p_z_position '+str(phantom_bead_1[2])+\
                            ' -var bead_2p_x_position '+str(phantom_bead_2[0])+' -var bead_2p_y_position '+str(phantom_bead_2[1])+' -var bead_2p_z_position '+str(phantom_bead_2[2])+' -var bead_3p_x_position '+str(phantom_bead_3[0])+\
                            ' -var bead_3p_y_position '+str(phantom_bead_3[1])+' -var bead_3p_z_position '+str(phantom_bead_3[2])+\
                            ' -var angle_stiff '+str(bending_stiffness[0])+' -var spring_stiffness '+str(internal_stiffness[m])+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+\
                            ' -var VP_ave_freq '+str(VP_ave_freq)+' -var realisation_index '+str(realisation_index_[j])+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+\
                            ' -var rand_int_2 '+rand_int_2+' -var rand_int_3 '+rand_int_3+' -var no_SRD '+no_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+\
                            ' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -in '+abs_path_2_lammps_script+' \n &'  #>> '+prod_run_file_name+' & \n'
                            
                            run_code_list.append(run_code_individual)
                            run_code=run_code +run_code_individual

                            
                            run_code = run_code[:-3]

                            py2bash_launch_overwriter.py2bash_launch_overwriter(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,num_proc,wd_path,extra_code,run_code,data_transfer_instructions)
    return run_code_list



# %%
run_code_list=sim_file_prod_flat_elastic_restarts_MYRIAD(restart_count,sim_batchcode,strain_start_point,restart_frequency,coordinates_tuple_3d,erate,equilibrium_triangle_side_length,var_choice_1,var_choice_2,internal_stiffness,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,realisation_index_,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,no_timestep_,thermo_freq,md_timestep,i_,j_,box_size_bar,box_size_index,Path_2_shell_scirpts,Path_2_generic,fluid_name)
    
# %%
#os.mkdir("/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/Simulation_run_folder/test_run_flat_elastic_100_realisations/")
os.chdir("/Volumes/Backup Plus 1/PhD_/Rouse Model simulations/Using LAMMPS imac/")
MyFile=open('initial_restart_test_'+str(strain_start_point)+'.sh','w')

for element in run_code_list:
     MyFile.write(element)
     MyFile.write('\n')
MyFile.close()

# %%
