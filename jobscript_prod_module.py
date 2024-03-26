
import os 
import datetime
import numpy as np 


def sim_file_prod_flat_elastic_initial_MYRIAD(restart_frequency,
                                              coordinates_tuple_3d,
                                              erate,
                                              equilibrium_triangle_side_length,
                                              var_choice_1,
                                              var_choice_2,
                                              internal_stiffness,
                                              bending_stiffness,
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
                                              i_,
                                              j_,
                                              box_size_bar,
                                              box_size_index,
                                              Path_2_shell_scirpts,
                                              Path_2_generic,
                                              fluid_name,
                                              total_strain,
                                              collision_time_negative_bar,
                                              srd_count):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_eqts_'+\
        str(equilibrium_triangle_side_length)+'_realisations_'+str(j_)+'_box_size_'+\
        str(box_size_bar[box_size_index])+'_bendstiff_'+str(bending_stiffness[0])+'_intstiff_'+\
        str(internal_stiffness[0])+'_'+str(internal_stiffness[-1])+'_'+ META_DATA
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
                            simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+\
                            '_realisation_'+str(j)+'_Bk_'+str(bending_stiffness[0])+'_np_'+str(np_req[k])+\
                            '_no_timesteps_'+str(no_timestep_[k])+'_intstiff_'+str(var_choice_2[m])+\
                            '_eqsl_'+str(equilibrium_triangle_side_length)+'_erate_'+str(erate[k])+'_'
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