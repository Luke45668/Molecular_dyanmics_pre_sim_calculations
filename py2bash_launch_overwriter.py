#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 14:57:00 2023

This script will be used to overwrite variables in my generic myriad launch script 
@author: lukedebono
"""
import regex as re
import os 
from datetime import datetime
#import shlex as sh 


#%% Read in the generic script 

#Path_2_shell=input('Please insert  absolute path to generic launch script')
def py2bash_launch_overwriter_mpi(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions):
    #Path_2_generic='/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD/' 
    os.chdir(Path_2_generic)
    filename_generic_launcher = 'generic_myriad_launch.sh'
    
    
    with open(filename_generic_launcher,'r+') as f:
    
           generic_myriad_launch_string = str(f.readlines())
           print(generic_myriad_launch_string)
           
    f.close()


    
    #%% use regex to replace key word and produce launch script 
    
    #change to batch folder
    os.chdir(simulation_batch_folder)
    
   
    #simulation_run_name = 'simple_LAMMPS_job' # always have a space at end
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    
    
    #JOBNAME 
    Specific_JOBNAME=simulation_run_name +META_DATA
    JOBNAME_regex = re.compile(r'JOBNAME')
    specific_myriad_launch_string = re.sub(JOBNAME_regex,Specific_JOBNAME,generic_myriad_launch_string) 
    
    #email 
    #specific_email = 'luke.debono.21@ucl.ac.uk'
    email_regex= re.compile(r'EMAIL')
    specific_myriad_launch_string = re.sub(email_regex,specific_email,specific_myriad_launch_string) 
    
    #priority phrase
    # swtich on to get gold back
    priority='' #'#$ -P Gold \n#$ -A hpc.13'
    priority_regex= re.compile(r'PRIORITY')
    specific_myriad_launch_string = re.sub(priority_regex,priority,specific_myriad_launch_string) 
    
    #wall clock
    #wall_time='00:02:00'
    wall_time_regex = re.compile(r'HH:MM:SS')
    specific_myriad_launch_string = re.sub(wall_time_regex,wall_time,specific_myriad_launch_string) 
    
    
    #ram request
    
    #ram_requirement='5G'
    ram_requirement_regex=re.compile(r'RAMMEM')
    specific_myriad_launch_string = re.sub(ram_requirement_regex,ram_requirement,specific_myriad_launch_string) 
    
    
    #tempdir space requested
    #tempdir_req='1G'
    tempdir_req_regex = re.compile(r'TMPMEM')
    specific_myriad_launch_string= re.sub(tempdir_req_regex, tempdir_req, specific_myriad_launch_string)
    
    #number of tasks 
    #num_task_req=''
    num_task_req_regex=re.compile(r'-t 1-TOTALNUMTASKS')
    specific_myriad_launch_string= re.sub(num_task_req_regex, num_task_req, specific_myriad_launch_string)
    
    
    
    #number of processors per task x the number of tasks 
    #='4'# could insert varaible from before 
    np_req_regex= re.compile(r'NTASKSPERNODE')
    specific_myriad_launch_string= re.sub(np_req_regex, np_req, specific_myriad_launch_string)
    #working directory path 
    #wd_path= '/home/ucahlrl/Scratch/output/' # clearly this needs to be more complicated 
    wd_path_regex= re.compile(r'LOCALDIR')
    specific_myriad_launch_string= re.sub(wd_path_regex, wd_path, specific_myriad_launch_string)
    
    #extra pre run code 
    #extra_code= 'module unload compilers mpi \n module load compilers/gnu/4.9.2 \n module load mpi/openmpi/1.10.1/gnu-4.9.2 \n module load cuda/7.5.18/gnu-4.9.2'
    extra_code_regex = re.compile(r'EXTRA')
    specific_myriad_launch_string= re.sub(extra_code_regex,extra_code , specific_myriad_launch_string)
    
    #run code 
    #run_code = 'mpirun -np 4 ./lmp_g++_openmpi -in in.srd.pure'
    run_code_regex= re.compile(r'RUN')
    specific_myriad_launch_string = re.sub(run_code_regex,run_code , specific_myriad_launch_string)
    
    #
    # transferring data to where ?
    data_transfer_instructions=''
    data_transfer_instructions_regex= re.compile(r'DATATRANSFER')
    specific_myriad_launch_string = re.sub(data_transfer_instructions_regex,data_transfer_instructions , specific_myriad_launch_string)
    specific_myriad_launch_string=specific_myriad_launch_string[1:]
    specific_myriad_launch_string=specific_myriad_launch_string[:-1]
    
    
    
    
    # get rid of quotation marks
    quote_mark='\'#'
    quote_mark_regex=re.compile(r'\ \'#')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
     
     
    
     # get rid of quotation marks
    quote_mark=''
    quote_mark_regex=re.compile(r'\'')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)

    quote_mark=''
    quote_mark_regex=re.compile(r',')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
    # this swap seems to make the script format itself properly 

    quote_mark='\\n'
    quote_mark_regex=re.compile(r'\\n')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
    
    
    
    #print(specific_myriad_launch_string)
    with  open(Specific_JOBNAME+'.sh','w') as f:
     
        
       
            # makes a new file for the job 
         f.write(specific_myriad_launch_string)
         f.close()
    #os.save(specific_myriad_launch_string,'.sh')

def py2bash_launch_overwriter_serial(Path_2_generic,simulation_batch_folder,simulation_run_name,specific_email,wall_time,ram_requirement,tempdir_req,num_task_req,np_req,wd_path,extra_code,run_code,data_transfer_instructions):
    #Path_2_generic='/Volumes/Backup Plus/PhD_/Rouse Model simulations/Using LAMMPS imac/Shell_scripts_for_MYRIAD/' 
    os.chdir(Path_2_generic)
    filename_generic_launcher = 'generic_myriad_launch_serial.sh'
    
    
    with open(filename_generic_launcher,'r+') as f:
    
           generic_myriad_launch_string = str(f.readlines())
           print(generic_myriad_launch_string)
           
    f.close()
    
    
    #%% use regex to replace key word and produce launch script 
    
    #change to batch folder
    os.chdir(simulation_batch_folder)
    
   
    #simulation_run_name = 'simple_LAMMPS_job' # always have a space at end
    META_DATA = str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))
    
    
    #JOBNAME 
    Specific_JOBNAME=simulation_run_name +META_DATA
    JOBNAME_regex = re.compile(r'JOBNAME')
    specific_myriad_launch_string = re.sub(JOBNAME_regex,Specific_JOBNAME,generic_myriad_launch_string) 
    
    #email 
    #specific_email = 'luke.debono.21@ucl.ac.uk'
    email_regex= re.compile(r'EMAIL')
    specific_myriad_launch_string = re.sub(email_regex,specific_email,specific_myriad_launch_string) 
    
    #priority phrase
    priority='#$ -P Gold \n#$ -A hpc.13'
    priority_regex= re.compile(r'PRIORITY')
    specific_myriad_launch_string = re.sub(priority_regex,priority,specific_myriad_launch_string) 
    
    #wall clock
    #wall_time='00:02:00'
    wall_time_regex = re.compile(r'HH:MM:SS')
    specific_myriad_launch_string = re.sub(wall_time_regex,wall_time,specific_myriad_launch_string) 
    
    
    #ram request
    
    #ram_requirement='5G'
    ram_requirement_regex=re.compile(r'RAMMEM')
    specific_myriad_launch_string = re.sub(ram_requirement_regex,ram_requirement,specific_myriad_launch_string) 
    
    
    #tempdir space requested
    #tempdir_req='1G'
    tempdir_req_regex = re.compile(r'TMPMEM')
    specific_myriad_launch_string= re.sub(tempdir_req_regex, tempdir_req, specific_myriad_launch_string)
    
    
    
    #number of processors per task x the number of tasks 
    #='4'# could insert varaible from before 
    np_req_regex= re.compile(r'NTASKSPERNODE')
    specific_myriad_launch_string= re.sub(np_req_regex, np_req, specific_myriad_launch_string)
    #working directory path 
    #wd_path= '/home/ucahlrl/Scratch/output/' # clearly this needs to be more complicated 
    wd_path_regex= re.compile(r'LOCALDIR')
    specific_myriad_launch_string= re.sub(wd_path_regex, wd_path, specific_myriad_launch_string)
    
    #extra pre run code 
    #extra_code= 'module unload compilers mpi \n module load compilers/gnu/4.9.2 \n module load mpi/openmpi/1.10.1/gnu-4.9.2 \n module load cuda/7.5.18/gnu-4.9.2'
    extra_code_regex = re.compile(r'EXTRA')
    specific_myriad_launch_string= re.sub(extra_code_regex,extra_code , specific_myriad_launch_string)
    
    #run code 
    #run_code = 'mpirun -np 4 ./lmp_g++_openmpi -in in.srd.pure'
    run_code_regex= re.compile(r'RUN')
    specific_myriad_launch_string = re.sub(run_code_regex,run_code , specific_myriad_launch_string)
    
    #
    # transferring data to where ?
    data_transfer_instructions=''
    data_transfer_instructions_regex= re.compile(r'DATATRANSFER')
    specific_myriad_launch_string = re.sub(data_transfer_instructions_regex,data_transfer_instructions , specific_myriad_launch_string)
    specific_myriad_launch_string=specific_myriad_launch_string[1:]
    specific_myriad_launch_string=specific_myriad_launch_string[:-1]
    
    
    
    
    # get rid of quotation marks
    quote_mark='\'#'
    quote_mark_regex=re.compile(r'\ \'#')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
     
     
    
     # get rid of quotation marks
    quote_mark=''
    quote_mark_regex=re.compile(r'\'')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)

    quote_mark=''
    quote_mark_regex=re.compile(r',')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
    # this swap seems to make the script format itself properly 

    quote_mark='\\n'
    quote_mark_regex=re.compile(r'\\n')
    specific_myriad_launch_string=re.sub(quote_mark_regex, quote_mark, specific_myriad_launch_string)
    
    
    
    #print(specific_myriad_launch_string)
    with  open(Specific_JOBNAME+'.sh','w') as f:
     
        
       
            # makes a new file for the job 
         f.write(specific_myriad_launch_string)
         f.close()
    #os.save(specific_myriad_launch_string,'.sh')