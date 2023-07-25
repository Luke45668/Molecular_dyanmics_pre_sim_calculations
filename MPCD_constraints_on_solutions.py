#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 12:12:41 2023
This function applies the constraints discussed in the MPCD chapter on mapping from a real fluid. 
This could be optimised by using indirect BOOLEAN if statements look up bookmark on google. 

@author: lukedebono
"""
import numpy as np 
def MPCD_constraints(no_timesteps,min_particle_count,sc_neg_soln,sc_pos_soln,srd_ratio_tolerance,max_particle_count,number_SRD_particles_wrt_pf_cp_mthd_1_pos,number_SRD_particles_wrt_pf_cp_mthd_1_neg,mean_free_path_pf_SRD_particles_cp_mthd_1_neg,mean_free_path_pf_SRD_particles_cp_mthd_1_pos,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos,Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg,Solvent_bead_SRD_box_density_cp_1,tolerance,SRD_box_size_wrt_solid_beads,comparison_pos,comparison_neg):
    
      
      
      
    for z in range(0,Solvent_bead_SRD_box_density_cp_1.size):
        for i in range(0,SRD_box_size_wrt_solid_beads.size):
    # for z in range(0,1):
    #     for i in  range(0,1):
            
        
            
            #### Need to add SRD/MD ratio constraint 
            # #positive 
            if Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]/no_timesteps > 0.1:
               
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] =float("NAN") 
                #fail_count_pos_constraint_total_time=fail_count_pos_constraint_total_time+1
            else:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] =Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]
            # negative 
            if Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i]/no_timesteps > 0.1:
                #print('Solution Fail! 1')
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] =float("NAN") 
                #fail_count_neg_constraint_total_time=fail_count_neg_constraint_total_time+1
            else:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] =Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i]
            
            #print(number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i])    
            ##### checking the solution isnt too far from integer value of SRD/MD ratio 
            
            if np.abs(comparison_pos[z,i])>tolerance:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] =float("NAN")    
                #fail_count_pos_constraint_solution_tolerance=fail_count_pos_constraint_solution_tolerance+1
            else:
             #print()
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]=Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]
               
                                
            if np.abs(comparison_neg[z,i])>tolerance:
                #print('Solution Fail 2')
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] =float("NAN")   
                ##fail_count_neg_constraint_solution_tolerance=fail_count_neg_constraint_solution_tolerance+1
            else:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i]=Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i]
            # print()
            #print(number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i])     
            #Knudsen number constraint
            
            #negative solution   
                        
            if (mean_free_path_pf_SRD_particles_cp_mthd_1_neg[z,i]/SRD_box_size_wrt_solid_beads[i])< 0.1:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] =Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] 
            else:
                     #mean_free_path_pf_SRD_particles_cp_mthd_1_neg[z,i]= float("NAN")
                #print('Solution Fail! 3')    
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] = float("NAN")
            
                number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i]=float("NAN")
                
                #fail_count_neg_constraint_Kn_number=fail_count_neg_constraint_Kn_number+1
                
                     #print("Knudsen number NOT in the continuum regime")
            #print(number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i])             
            #positive solution 
            
            if (mean_free_path_pf_SRD_particles_cp_mthd_1_pos[z,i]/SRD_box_size_wrt_solid_beads[i])< 0.1:
                     #print("Knudsen number in the continuum regime")
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]=Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]
            else:
                   
                     Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] = float("NAN")
                
                     #number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i]=float("NAN")
                
                #fail_count_pos_constraint_Kn_number=fail_count_pos_constraint_Kn_number+1
                
                     #print("Knudsen number NOT in the continuum regime")
            #print(number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i])                    
            # testing for SRd particle count max 
            
            # positive       
                     
            if number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i] <max_particle_count:
               number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i]=number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i]
            else: 
              
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] = float("NAN")
                #print("All tests with SRD count > "+str(max_particle_count)+" removed ")
               
                mean_free_path_pf_SRD_particles_cp_mthd_1_pos[z,i]= float("NAN")
                #fail_count_pos_max_particle_count=fail_count_pos_max_particle_count+1
                
            # negative 
            #print(number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i])    
            if number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i] <max_particle_count:
                
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] =Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i]     
            else: 
                
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] = float("NAN")
                
                #print('Solution Fail!4')
               
               # print("All tests with SRD count > "+str(max_particle_count)+" removed ")
                mean_free_path_pf_SRD_particles_cp_mthd_1_neg[z,i]= float("NAN")
                #fail_count_neg_max_particle_count=fail_count_neg_max_particle_count+1
                
                
            # testing for min particle count 
            
            #positive 
                
            if number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i] >min_particle_count:
                number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i]=number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i]   
            else: 
              
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] = float("NAN")
                #print("All tests with SRD count < "+str(min_particle_count)+" removed ")
               
                mean_free_path_pf_SRD_particles_cp_mthd_1_pos[z,i]= float("NAN")
                #fail_count_pos_min_particle_count=fail_count_pos_min_particle_count+1
                
            
            #negative     
                
            if number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i] >min_particle_count:
                
                number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i]=number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i]  
            else: 
                #print('Solution Fail!5')
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] = float("NAN")
               
                #print("All tests with SRD count < "+str(min_particle_count)+" removed ")
                mean_free_path_pf_SRD_particles_cp_mthd_1_neg[z,i]= float("NAN")
                #fail_count_neg_min_particle_count=fail_count_neg_min_particle_count+1
                
                
            ## SRD ratio tolerance test
            
            #negative 
                 
            if Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] <srd_ratio_tolerance:
                #print('Solution Fail!6')
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] = float("NAN")
                
                number_SRD_particles_wrt_pf_cp_mthd_1_neg[z,i]=float("NAN")
                #print("SRD MD Ratio <"+str(srd_ratio_tolerance))    
                #fail_count_neg_SRD_ratio_tolerance=fail_count_neg_SRD_ratio_tolerance+1
            else: 
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i]=Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i]
            
            #positive     
                
            if Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] <srd_ratio_tolerance:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] = float("NAN")
                
                number_SRD_particles_wrt_pf_cp_mthd_1_pos[z,i]=float("NAN")
                #print("SRD MD Ratio <"+str(srd_ratio_tolerance))    
                #fail_count_pos_SRD_ratio_tolerance=fail_count_pos_SRD_ratio_tolerance+1
            else: 
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] = Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] 
                
            #Schmidt number test 
            
            #positive 
            
            if sc_neg_soln[z,i]<100: # was 100 but wagner paper says 13 is good 
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] = float("NAN")
                mean_free_path_pf_SRD_particles_cp_mthd_1_neg[z,i]= float("NAN")
                #fail_count_pos_Sc_num=fail_count_pos_Sc_num+1
            else:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] =Number_MD_steps_per_SRD_with_pf_cp_mthd_1_neg[z,i] 
            
            # negative 
            
            if sc_pos_soln[z,i]<100: # was 100 but wagner pape says 13 is good 
                 #print('Solution Fail!7')
                 Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i] = float("NAN")
                 mean_free_path_pf_SRD_particles_cp_mthd_1_pos[z,i]= float("NAN")
                 #fail_count_neg_Sc_num=fail_count_neg_Sc_num+1
            else:
                Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]= Number_MD_steps_per_SRD_with_pf_cp_mthd_1_pos[z,i]
            
      