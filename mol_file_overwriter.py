#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/03/24 

This script will be used to overwrite generic molecule files
@author: lukedebono
"""
import regex as re
import os 
from datetime import datetime
#import shlex as sh 


#%% Read in the generic script 

#Path_2_shell=input('Please insert  absolute path to generic launch script')
def mol_file_coordinate_overwriter(Path_2_generic,filename_generic,simulation_batch_folder,coord_string,name_particle,rounded_normal):
    #Path_2_generic='/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD/' 
    os.chdir(Path_2_generic)
    #filename_generic_launcher = 'generic_myriad_launch.sh'
    
    
    with open(filename_generic,'r+') as f:
    
           generic_mol_file_string = str(f.readlines())
           print(generic_mol_file_string)
           
    f.close()
    
    
    # use regex to replace key word and produce launch script 
    
    #change to batch folder
    os.chdir(simulation_batch_folder)

    
    #INSERT_COORDS
    Specific_coords=coord_string 
    coords_regex = re.compile(r'INSERT_COORDS')
    specific_myriad_launch_string = re.sub(coords_regex,Specific_coords,generic_mol_file_string) 
    

    
    
    #get rid of quotation marks
    quote_mark='\'#'
    quote_mark_regex=re.compile(r'\ \'#')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
     
     
    
#      #get rid of quotation marks
    quote_mark=''
    quote_mark_regex=re.compile(r'\'')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)

    quote_mark=''
    quote_mark_regex=re.compile(r',')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
    
    # remove leading white space
    quote_mark=''
    quote_mark_regex=re.compile(r'^[\s]+')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
    # this swap seems to make the script format itself properly 
    quote_mark='\\n'
    quote_mark_regex=re.compile(r'\\n')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
    

    specific_myriad_launch_string=specific_myriad_launch_string[1:]
    specific_myriad_launch_string=specific_myriad_launch_string[:-1]
    #specific_myriad_launch_string=specific_myriad_launch_string.lstrip(' ')
   # print(specific_myriad_launch_string)
    Specific_molfile_name=name_particle+'_'+str(rounded_normal[0])+'_'+str(rounded_normal[1])+'_'+str(rounded_normal[1])
    
    #print(specific_myriad_launch_string)
    with  open('mol.'+Specific_molfile_name,'w') as f:
     
        
       
            # makes a new file for the job 
         f.write(specific_myriad_launch_string)
         f.close()
    #os.save(specific_myriad_launch_string,'.sh')
         
    return specific_myriad_launch_string




