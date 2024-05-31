#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:26:55 2023
Function to check the box size constraint is satisified and the simulation box will have an integer number of bins
@author: lukedebono
"""
def box_size_dim_and_integer_SRD_bin_count_constraint_func(k_b,T_K,gamma_1,rho_s,nu_s,box_side_length_scaled,SRD_box_size_wrt_solid_beads,SRD_box_size_wrt_solid_beads_check):
    Box_size_dimensional_constraint = (k_b*T_K*gamma_1)/18*rho_s*nu_s**(2)
    
    if SRD_box_size_wrt_solid_beads_check.any() < Box_size_dimensional_constraint.any():#is this correct 
        print("Constraint 1 on box size NOT satisfied")
    else:
        print("Constraint 1 on box size satisfied")
     
    number_of_SRD_bins_per_dimension= box_side_length_scaled/SRD_box_size_wrt_solid_beads
    for i in range(0,SRD_box_size_wrt_solid_beads.size):
     if number_of_SRD_bins_per_dimension[0,i].is_integer():
        print("Integer number of bins achieved")
     else:
        print("Non-integer number of bins will throw error")
        
    return Box_size_dimensional_constraint
