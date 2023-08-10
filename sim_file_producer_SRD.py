#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:49:47 2023

This is my function which takes acceptable solutions from the SRD calc and produces 
a file named specifically for the the simulation run to store the bash scripts. 




@author: lukedebono
"""

import numpy as np
import os 
from datetime import datetime
import py2bash_launch_overwriter
import py2bash_launch_overwriter_kathleen

# note on how to use this function: just change the swap rate and swap number vector , need to add a check to test whether the correct amount of cores has been accquired. 
def sim_file_prod_neg_soln(phi,solution_choice_tuple,lengthscale_parameter,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,prod_run_file_name,realisation_index_,equilibration_timesteps,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,num_proc,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,scaled_timestep,mass_fluid_particle_wrt_pf_cp_mthd_1,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_neg,swap_number,i_,j_,swap_rate,box_side_length_scaled,scaled_temp,eta_s,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_phi_'+phi+'_validation_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_no_timesteps_'+str(no_timesteps)+'_tstep_'+str(scaled_timestep)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    # test to check consistency of cores request 
    
    
    for j in range(i_,j_): #or now just use one realisation 
    
        for k in range(0,swap_number.size):
                param_set_code=str(np.random.randint(0, 1000000))
                #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_''_swap_number_'+str(swap_number[k])+'_test_run_'
                simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                run_code=''   
                for m in range(0,swap_rate.size):#range(0,1):  
                    #or l in range(0,array_entry.size):
                        #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                    
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1_neg)) 
                        #print(no_SRD)
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1)
                        box_size = str(box_side_length_scaled)
                        timestep_input= str(scaled_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)
                        SRD_MD_ratio=str(int(SRD_MD_ratio))
                        lamda= str(mean_free_path_pf_SRD_particles_cp_mthd_1_neg )
                        grid_size = str(SRD_box_size_wrt_solid_beads)
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps = str(no_timesteps)
                        temp_=str(scaled_temp)
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))


                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+'  -in '+abs_path_2_lammps_script+' & \n'  #>> '+prod_run_file_name+' & \n'
                        run_code=run_code +run_code_individual

                run_code = run_code[:-3]

                py2bash_launch_overwriter.py2bash_launch_overwriter(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)


def sim_file_prod_neg_soln_individual(phi,solution_choice_tuple,lengthscale_parameter,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,prod_run_file_name,realisation_index_,equilibration_timesteps,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,num_proc,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,scaled_timestep,mass_fluid_particle_wrt_pf_cp_mthd_1,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_neg,swap_number,i_,j_,swap_rate,box_side_length_scaled,scaled_temp,eta_s,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_phi_'+phi+'_validation_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_no_timesteps_'+str(no_timesteps)+'_tstep_'+str(scaled_timestep)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    # test to check consistency of cores request 
    
    
    for j in range(i_,j_): #or now just use one realisation 
    
        for k in range(0,swap_number.size):
                param_set_code=str(np.random.randint(0, 1000000))
                #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_''_swap_number_'+str(swap_number[k])+'_test_run_'
                # simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                # run_code=''   
                for m in range(0,swap_rate.size):#range(0,1):  
                    #or l in range(0,array_entry.size):
                        #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                        simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                        run_code=''
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1_neg)) 
                        #print(no_SRD)
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1)
                        box_size = str(box_side_length_scaled)
                        timestep_input= str(scaled_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)
                        SRD_MD_ratio=str(int(SRD_MD_ratio))
                        lamda= str(mean_free_path_pf_SRD_particles_cp_mthd_1_neg )
                        grid_size = str(SRD_box_size_wrt_solid_beads)
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps = str(no_timesteps)
                        temp_=str(scaled_temp)
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))


                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+'  -in '+abs_path_2_lammps_script+' & \n'  #>> '+prod_run_file_name+' & \n'
                        run_code=run_code +run_code_individual

                        run_code = run_code[:-3]

                        py2bash_launch_overwriter.py2bash_launch_overwriter(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)

def sim_file_prod_neg_soln_solid_inc_individual(spring_constant,phi,mass_solid,particle_x_upper_nd,particle_y_upper_nd,particle_z_upper_nd,particle_x_lower_nd,particle_y_lower_nd,particle_z_lower_nd,solution_choice_tuple,lengthscale_parameter,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,prod_run_file_name,realisation_index_,equilibration_timesteps,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,num_proc,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,scaled_timestep,mass_fluid_particle_wrt_pf_cp_mthd_1,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_neg,swap_number,i_,j_,swap_rate,box_side_length_scaled,scaled_temp,eta_s,Path_2_shell_scirpts,Path_2_generic,fluid_name,r_particle_scaled):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_phi_'+phi+'_validation_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_no_timesteps_'+str(no_timesteps)+'_tstep_'+str(scaled_timestep)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    # test to check consistency of cores request 
    
    
    for j in range(i_,j_): #or now just use one realisation 
    
        for k in range(0,swap_number.size):
                param_set_code=str(np.random.randint(0, 1000000))
                #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_''_swap_number_'+str(swap_number[k])+'_test_run_'
                # simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                # run_code=''   
                for m in range(0,swap_rate.size):#range(0,1):  
                    #or l in range(0,array_entry.size):
                        #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                        
                        simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_k_'+str(spring_constant[m])+'_'
                        run_code=''
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1_neg)) 
                        #print(no_SRD)
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1)
                        box_size = str(box_side_length_scaled)
                        timestep_input= str(scaled_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)
                        SRD_MD_ratio=str(int(SRD_MD_ratio))
                        lamda= str(mean_free_path_pf_SRD_particles_cp_mthd_1_neg )
                        grid_size = str(SRD_box_size_wrt_solid_beads)
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps = str(no_timesteps)
                        temp_=str(scaled_temp)
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))
                        


                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var spring_constant '+str(spring_constant[m])+' -var fluid_name '+fluid_name +' -var mass_solid '+mass_solid+' -var particle_x_upper_nd '+particle_x_upper_nd+' -var particle_y_upper_nd '+particle_y_upper_nd+' -var particle_z_upper_nd '+particle_z_upper_nd +' -var  particle_x_lower_nd '+particle_x_lower_nd+' -var particle_y_lower_nd '+particle_y_lower_nd+' -var particle_z_lower_nd '+particle_z_lower_nd+' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -var r_particle '+str(r_particle_scaled)+' -in '+abs_path_2_lammps_script+' & \n'  #>> '+prod_run_file_name+' & \n'
                        run_code=run_code +run_code_individual


                        run_code = run_code[:-3]

                        py2bash_launch_overwriter.py2bash_launch_overwriter(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)



def sim_file_prod_neg_soln_solid_inc(phi,mass_solid,particle_x_upper_nd,particle_y_upper_nd,particle_z_upper_nd,particle_x_lower_nd,particle_y_lower_nd,particle_z_lower_nd,solution_choice_tuple,lengthscale_parameter,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,prod_run_file_name,realisation_index_,equilibration_timesteps,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,num_proc,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,scaled_timestep,mass_fluid_particle_wrt_pf_cp_mthd_1,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_neg,swap_number,i_,j_,swap_rate,box_side_length_scaled,scaled_temp,eta_s,Path_2_shell_scirpts,Path_2_generic,fluid_name,r_particle_scaled):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_phi_'+phi+'_validation_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_no_timesteps_'+str(no_timesteps)+'_tstep_'+str(scaled_timestep)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    
    
    for j in range(i_,j_): #or now just use one realisation 
    
        for k in range(0,swap_number.size):
                param_set_code=str(np.random.randint(0, 1000000))
                #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_''_swap_number_'+str(swap_number[k])+'_test_run_'
                simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                run_code=''   
                for m in range(0,swap_rate.size):#range(0,1):  
                    #or l in range(0,array_entry.size):
                        #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                    
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1_neg)) 
                        #print(no_SRD)
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1)
                        box_size = str(box_side_length_scaled)
                        timestep_input= str(scaled_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)
                        SRD_MD_ratio=str(int(SRD_MD_ratio))
                        lamda= str(mean_free_path_pf_SRD_particles_cp_mthd_1_neg )
                        grid_size = str(SRD_box_size_wrt_solid_beads)
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps = str(no_timesteps)
                        temp_=str(scaled_temp)
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))

                        #particle_x_upper_nd,particle_y_upper_nd,particle_z_upper_nd,particle_x_lower_nd,particle_y_lower_nd,particle_z_lower_nd

                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var fluid_name '+fluid_name +' -var mass_solid '+mass_solid+' -var particle_x_upper_nd '+particle_x_upper_nd+' -var particle_y_upper_nd '+particle_y_upper_nd+' -var particle_z_upper_nd '+particle_z_upper_nd +' -var  particle_x_lower_nd '+particle_x_lower_nd+' -var particle_y_lower_nd '+particle_y_lower_nd+' -var particle_z_lower_nd '+particle_z_lower_nd+' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -var r_particle '+str(r_particle_scaled)+' -in '+abs_path_2_lammps_script+' & \n'  #>> '+prod_run_file_name+' & \n'
                        run_code=run_code +run_code_individual

                run_code = run_code[:-3]

                py2bash_launch_overwriter.py2bash_launch_overwriter(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)


    
def sim_file_prod_neg_soln_kathleen(phi,hypthread,solution_choice_tuple,lengthscale_parameter,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,prod_run_file_name,realisation_index_,equilibration_timesteps,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,num_proc,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,scaled_timestep,mass_fluid_particle_wrt_pf_cp_mthd_1,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_neg,swap_number,i_,j_,swap_rate,box_side_length_scaled,scaled_temp,eta_s,Path_2_shell_scirpts,Path_2_generic,fluid_name):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_phi_'+phi+'_validation_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_no_timesteps_'+str(no_timesteps)+'_tstep_'+str(scaled_timestep)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    
    
    for j in range(i_,j_): #or now just use one realisation 
    
        for k in range(0,swap_number.size):
                param_set_code=str(np.random.randint(0, 1000000))
                #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_''_swap_number_'+str(swap_number[k])+'_test_run_'
                simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                run_code=''   
                for m in range(0,swap_rate.size):#range(0,1):  
                    #or l in range(0,array_entry.size):
                        #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                    
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1_neg)) 
                        #print(no_SRD)
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1)
                        box_size = str(box_side_length_scaled)
                        timestep_input= str(scaled_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)
                        SRD_MD_ratio=str(int(SRD_MD_ratio))
                        lamda= str(mean_free_path_pf_SRD_particles_cp_mthd_1_neg )
                        grid_size = str(SRD_box_size_wrt_solid_beads)
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps = str(no_timesteps)
                        temp_=str(scaled_temp)
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))


                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var fluid_name '+fluid_name +' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+'  -in '+abs_path_2_lammps_script+' & \n'  #>> '+prod_run_file_name+' & \n'
                        run_code=run_code +run_code_individual

                run_code = run_code[:-3]

                py2bash_launch_overwriter_kathleen.py2bash_launch_overwriter(hypthread,Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)


def sim_file_prod_neg_soln_solid_inc_kathleen(phi,hypthread,mass_solid,particle_x_upper_nd,particle_y_upper_nd,particle_z_upper_nd,particle_x_lower_nd,particle_y_lower_nd,particle_z_lower_nd,solution_choice_tuple,lengthscale_parameter,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,prod_run_file_name,realisation_index_,equilibration_timesteps,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,num_proc,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,scaled_timestep,mass_fluid_particle_wrt_pf_cp_mthd_1,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_neg,swap_number,i_,j_,swap_rate,box_side_length_scaled,scaled_temp,eta_s,Path_2_shell_scirpts,Path_2_generic,fluid_name,r_particle_scaled):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_phi_'+phi+'_validation_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_no_timesteps_'+str(no_timesteps)+'_tstep_'+str(scaled_timestep)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    
    
    for j in range(i_,j_): #or now just use one realisation 
    
        for k in range(0,swap_number.size):
                param_set_code=str(np.random.randint(0, 1000000))
                #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_''_swap_number_'+str(swap_number[k])+'_test_run_'
                simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                run_code=''   
                for m in range(0,swap_rate.size):#range(0,1):  
                    #or l in range(0,array_entry.size):
                        #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                    
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1_neg)) 
                        #print(no_SRD)
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1)
                        box_size = str(box_side_length_scaled)
                        timestep_input= str(scaled_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)
                        SRD_MD_ratio=str(int(SRD_MD_ratio))
                        lamda= str(mean_free_path_pf_SRD_particles_cp_mthd_1_neg )
                        grid_size = str(SRD_box_size_wrt_solid_beads)
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps = str(no_timesteps)
                        temp_=str(scaled_temp)
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))

                        #particle_x_upper_nd,particle_y_upper_nd,particle_z_upper_nd,particle_x_lower_nd,particle_y_lower_nd,particle_z_lower_nd

                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var fluid_name '+fluid_name +' -var mass_solid '+mass_solid+' -var particle_x_upper_nd '+particle_x_upper_nd+' -var particle_y_upper_nd '+particle_y_upper_nd+' -var particle_z_upper_nd '+particle_z_upper_nd +' -var  particle_x_lower_nd '+particle_x_lower_nd+' -var particle_y_lower_nd '+particle_y_lower_nd+' -var particle_z_lower_nd '+particle_z_lower_nd+' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -var r_particle '+str(r_particle_scaled)+' -in '+abs_path_2_lammps_script+' & \n'  #>> '+prod_run_file_name+' & \n'
                        run_code=run_code +run_code_individual

                run_code = run_code[:-3]

                py2bash_launch_overwriter_kathleen.py2bash_launch_overwriter(hypthread,Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)

def sim_file_prod_neg_soln_solid_inc_individual_kathleen(spring_constant,phi,hypthread,mass_solid,particle_x_upper_nd,particle_y_upper_nd,particle_z_upper_nd,particle_x_lower_nd,particle_y_lower_nd,particle_z_lower_nd,solution_choice_tuple,lengthscale_parameter,data_transfer_instructions,extra_code,wd_path,np_req,num_task_req,tempdir_req,wall_time,ram_requirement,prod_run_file_name,realisation_index_,equilibration_timesteps,VP_ave_freq,abs_path_2_lammps_exec,abs_path_2_lammps_script,num_proc,no_timesteps,thermo_freq,dump_freq,SRD_box_size_wrt_solid_beads,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,scaled_timestep,mass_fluid_particle_wrt_pf_cp_mthd_1,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,number_SRD_particles_wrt_pf_cp_mthd_1_neg,swap_number,i_,j_,swap_rate,box_side_length_scaled,scaled_temp,eta_s,Path_2_shell_scirpts,Path_2_generic,fluid_name,r_particle_scaled):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_phi_'+phi+'_validation_fluid_visc_'+str(eta_s)+'_temp_'+str(scaled_temp)+'_box_size_'+str(box_side_length_scaled)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_no_timesteps_'+str(no_timesteps)+'_tstep_'+str(scaled_timestep)+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    # test to check consistency of cores request 
    
    
    for j in range(i_,j_): #or now just use one realisation 
    
        for k in range(0,swap_number.size):
               
                #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_''_swap_number_'+str(swap_number[k])+'_test_run_'
                # simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_range_'+str(swap_rate[0])+'_'+str(swap_rate[swap_rate.size-1])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                # run_code=''   
                for m in range(0,swap_rate.size):#range(0,1):  
                    #or l in range(0,array_entry.size):
                        #simulation_run_name=fluid_name+'_val_param_sweep_solution_'+str(z)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                        param_set_code=str(np.random.randint(0, 1000000))
                        simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+'_'+str(lengthscale_parameter)+'_val_param_sweep_solution_'+str(solution_choice_tuple)+'_realisation_'+str(j)+'_swap_rate_'+str(swap_rate[m])+'_swap_number_'+str(swap_number[k])+'_test_run_'
                        run_code=''
                        no_SRD=str(int(number_SRD_particles_wrt_pf_cp_mthd_1_neg)) 
                        #print(no_SRD)
                        mass_SRD =str(mass_fluid_particle_wrt_pf_cp_mthd_1)
                        box_size = str(box_side_length_scaled)
                        timestep_input= str(scaled_timestep)
                        chunk = 20 # number of chunks to use for VP averaging
                        SRD_MD_ratio=np.round(Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg)
                        SRD_MD_ratio=str(int(SRD_MD_ratio))
                        lamda= str(mean_free_path_pf_SRD_particles_cp_mthd_1_neg )
                        grid_size = str(SRD_box_size_wrt_solid_beads)
                        dump_freq=str(dump_freq) # if you change the timestep rememebr to chaneg this 
                        thermo_freq = str(thermo_freq) # if you change the timestep rememebr to chaneg this 
                        no_timesteps = str(no_timesteps)
                        temp_=str(scaled_temp)
                        rand_int =str(np.random.randint(0, 1000000))
                        rand_int_1 =str( np.random.randint(0, 1000000))
                        rand_int_2 =str(np.random.randint(0, 1000000))


                        run_code_individual ='mpirun -np '+str(num_proc)+'  '+abs_path_2_lammps_exec+' -var spring_constant '+str(spring_constant[m])+' -var fluid_name '+fluid_name +' -var mass_solid '+mass_solid+' -var particle_x_upper_nd '+particle_x_upper_nd+' -var particle_y_upper_nd '+particle_y_upper_nd+' -var particle_z_upper_nd '+particle_z_upper_nd +' -var  particle_x_lower_nd '+particle_x_lower_nd+' -var particle_y_lower_nd '+particle_y_lower_nd+' -var particle_z_lower_nd '+particle_z_lower_nd+' -var  sim_batchcode '+str(sim_batchcode)+' -var swap_rate '+str(swap_rate[m])+' -var swap_number '+str(swap_number[k])+' -var VP_ave_freq '+str(VP_ave_freq)+' -var equilibration_timesteps '+str(equilibration_timesteps)+' -var chunk '+str(chunk)+' -var grid_size '+grid_size+' -var realisation_index '+str(realisation_index_[j])+' -var temp_ '+temp_+' -var lambda '+str(lamda)+' -var rand_int '+rand_int+' -var rand_int_1 '+rand_int_1+' -var rand_int_2 '+rand_int_2+' -var no_SRD '+no_SRD+' -var mass_SRD '+mass_SRD+' -var box_size '+box_size+' -var timestep_input '+timestep_input+' -var SRD_MD_ratio '+SRD_MD_ratio+' -var dump_freq '+dump_freq+' -var thermo_freq '+thermo_freq+' -var no_timesteps '+no_timesteps+' -var r_particle '+str(r_particle_scaled)+' -in '+abs_path_2_lammps_script+' & \n'  #>> '+prod_run_file_name+' & \n'
                        run_code=run_code +run_code_individual


                        run_code = run_code[:-3]

                        py2bash_launch_overwriter_kathleen.py2bash_launch_overwriter(hypthread,Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions)

