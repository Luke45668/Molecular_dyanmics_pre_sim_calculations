##!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This code will produce the coordinates for an N sided elastic plate
"""
#%%
import os 
import numpy as np
import glob as glob
import matplotlib.pyplot as plt



#%% define function to calculate coords in 2d
original_triangle_side=3
rad_circle=np.sqrt(original_triangle_side)
n_side=5
internal_angle=(n_side-2)*np.pi 
internal_alpha=internal_angle/(2*n_side)
box_size=100
x=box_size/2
y=x

#set x and y to centre of plate 
def initial_coords_one_side(x,y,rad_circle,internal_alpha):
    stokes_1=np.array([x+(rad_circle*np.sin(internal_alpha)),
                       y-(rad_circle*np.cos(internal_alpha))])
    phantom_1=np.array([x,
                        y-(rad_circle*np.cos(internal_alpha))])
    stokes_2=np.array([x-(rad_circle*np.sin(internal_alpha)),
                       y-(rad_circle*np.cos(internal_alpha))])
    side_coords=np.array([stokes_1,
                          phantom_1,
                          stokes_2])
    return side_coords



def rotate_side(n_side,side_coords):


    theta=(2*np.pi/n_side)
    clockw_rotation=np.array([[np.cos(theta),-np.sin(theta)],
                              [np.sin(theta),np.cos(theta)]])
    return np.matmul(clockw_rotation,side_coords.T)


side_coords=initial_coords_one_side(x,y,
                                    rad_circle,
                                    internal_alpha)
rotated_coords=np.zeros((n_side,2,3))


i_side=1
input_coords=side_coords
rotated_coords[0,:,:]=input_coords.T
for i in range(1,5):
    

    rotated_coords[i,:,:]=rotate_side(n_side,input_coords)

    plt.plot(rotated_coords[i,0,:],rotated_coords[i,1,:],marker='x')
    i_side+=1
    input_coords=rotated_coords[i,:,:].T
plt.show()






    









#%% check shape 



#%% create molecule file 

