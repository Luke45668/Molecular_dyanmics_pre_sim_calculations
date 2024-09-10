import numpy as np
import os 
import datetime
from sim_file_producer_SRD import *
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

def out_put_freq_calc(no_steps,no_outputs):
    output_freq= no_steps/no_outputs
    if isinstance(output_freq, int) ==True :
            print("correct input")
    else:
            output_freq=np.floor(output_freq).astype('int')
           

    return output_freq

def compute_timesteps_for_strain(total_strain,erate,md_timestep,timestep_multiplier):
      no_timestep_=(np.round((total_strain/(erate*md_timestep*timestep_multiplier)),decimals=-3)).astype('int')

      return no_timestep_

def product_constraint_inputs(damp_upper_bound,damp_lower_bound,internal_stiff,initial_damp,n_sets):
      kdamp=internal_stiff*initial_damp
      Possible_damps=np.linspace(damp_lower_bound,damp_upper_bound,n_sets)

      new_internal_stiffness=kdamp/Possible_damps

      return new_internal_stiffness,Possible_damps




def sim_file_prod_flat_elastic_MYRIAD_all_erate_one_file(SRD_MD_ratio_,collision_time_negative_bar,bending_stiffness,damp,input_temp,
                                                         erate,
                                                         equilibrium_triangle_side_length,
                                                         var_choice_1,var_choice_2,internal_stiffness,
                                                         data_transfer_instructions,extra_code,wd_path,
                                                         np_req,num_task_req,tempdir_req,wall_time,
                                                         ram_requirement,realisation_index_,VP_ave_freq,
                                                         abs_path_2_lammps_exec,abs_path_2_lammps_script,
                                                         no_timestep_,thermo_freq,dump_freq,md_timestep,
                                                         i_,j_,box_size_bar,Path_2_shell_scirpts,Path_2_generic,fluid_name,timestep_multiplier,thermal_damp_multiplier):
    
    os.chdir(Path_2_shell_scirpts)
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    specific_email = 'luke.debono.21@ucl.ac.uk'
    simulation_batch_folder= 'simulation_batch_scripts_'+fluid_name+'_eqts_'+\
        str(equilibrium_triangle_side_length)+'_realisations_'+str(j_)+'_box_size_'+\
            str(box_size_bar)+'_bendstiff_'+str(bending_stiffness[0])+\
                '_intstiff_'+str(internal_stiffness[0])+'_'+str(internal_stiffness[-1])+\
                    '_tdamprange_'+str(thermal_damp_multiplier[0])+'_'+str(thermal_damp_multiplier[-1])+\
                    '_erate_'+str(var_choice_1[0])+'_'+str(var_choice_1[-1])+'_'+META_DATA
    os.mkdir(simulation_batch_folder)
    sim_batchcode=str(np.random.randint(0, 1000000))
    run_code_list=[]
    # test to check consistency of cores request 

    

    #for n in range(0,np_req.size):
        #or now just use one realisation 
    for h in range(thermal_damp_multiplier.size):
       thermal_damp_multi=thermal_damp_multiplier[h]

       for m in range(0,var_choice_2.size): 
             for k in range(erate.size):    
        #for k in range(0,1):   
            
                        #for m in range(0,1):  
                                for j in range(i_,j_):
                                    param_set_code=str(np.random.randint(0, 1000000))
                                    simulation_run_name=fluid_name+'_'+str(sim_batchcode)+'_'+param_set_code+\
                                        '_realisation_'+str(j)+'_Bk_'+str(bending_stiffness[0])+'_np_'+str(np_req)+\
                                            '_no_timesteps_'+str(no_timestep_[m,k])+'_intstiff_'+str(var_choice_2[m])+\
                                                '_eqsl_'+str(equilibrium_triangle_side_length)+'_erate_'+\
                                                    str(var_choice_1[k])+'_tdamp_'+str(thermal_damp_multi)+'_'
                                    run_code=''
                                
                                    #print(no_SRD)
                                    box_size = str(box_size_bar)
                                    timestep_input= str(md_timestep)
                                    SRD_MD_ratio=str(int(SRD_MD_ratio_))
                                    lamda= str(collision_time_negative_bar)
                                    dump_freq_=str(dump_freq[m,k])
                                    thermo_freq_ = str(thermo_freq[m,k])
                                    no_timesteps = str(no_timestep_[m,k])
                                    rand_int =str(np.random.randint(0, 1000000))
                                    rand_int_1 =str( np.random.randint(0, 1000000))
                                    rand_int_2 =str(np.random.randint(0, 1000000))
                                    rand_int_3=str(np.random.randint(0,1000000))
                                    
                                

                                                
                                    run_code_individual ="mpirun -np "+np_req+" "+abs_path_2_lammps_exec+' -var tdamp_multi '+str(thermal_damp_multi)+' -var timestep_multiplier '+str(timestep_multiplier[m,k])+' -var temp '+str(input_temp[k])+' -var damp '\
                                        +str(damp[0])+' -var erate_in '+str(erate[k])+' -var equilirbium_triangle_side_length '\
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
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_imac_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps-23Jun2022_with_SRD_pol/build_macbook_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_macbook_h5md/lmp_mpi'
#abs_path_2_lammps_exec='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/lammps_hirotori/build_m3maxbook/lmp_mpi'
#abs_path_2_lammps_script='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/LSC/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5'
#abs_path_2_lammps_script='/Volumes/Backup\ Plus\ 1/PhD_/Rouse\ Model\ simulations/Using\ LAMMPS\ imac/LSC/in.MPCD_with_hookean_flat_elastic_particle_only_dump_hdf5_chain'

