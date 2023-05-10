#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 09:23:00 2023
This function calculates the MPCD-SRD(-a) gamma functions 
@author: lukedebono
"""
import numpy as np
def SRD_a_minus_gamma_funcs_(Solvent_bead_SRD_box_density_cp_1):
    
    gamma_1 = ((1 - ((1-np.exp(-(Solvent_bead_SRD_box_density_cp_1)))/(Solvent_bead_SRD_box_density_cp_1))).T)
    gamma_2 = ( (((Solvent_bead_SRD_box_density_cp_1)+2)/((Solvent_bead_SRD_box_density_cp_1)-1)).T)
    return gamma_1,gamma_2
