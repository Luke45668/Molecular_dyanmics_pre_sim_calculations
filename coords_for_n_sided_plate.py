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
# n side should be odd
n_side=3
internal_angle=(n_side-2)*np.pi 
internal_alpha=internal_angle/(2*n_side)
box_size=100
box_centre=box_size/2

#set x and y to centre of plate 
def initial_coords_one_side(rad_circle,internal_alpha):
    stokes_1=np.array([(rad_circle*np.cos(internal_alpha)),
                       -(rad_circle*np.sin(internal_alpha))])
    phantom_1=np.array([0,
                        -(rad_circle*np.sin(internal_alpha))])
    stokes_2=np.array([-(rad_circle*np.cos(internal_alpha)),
                       -(rad_circle*np.sin(internal_alpha))])
    side_coords=np.array([stokes_1,
                          phantom_1,
                          stokes_2])
    return side_coords



def rotate_side(n_side,side_coords):


    theta=(2*np.pi/n_side)
    clockw_rotation=np.array([[np.cos(theta),-np.sin(theta)],
                              [np.sin(theta),np.cos(theta)]])
    return np.matmul(clockw_rotation,side_coords.T)


side_coords=initial_coords_one_side(
                                    rad_circle,
                                    internal_alpha)
rotated_coords=np.zeros((n_side,2,3))


i_side=1
input_coords=side_coords
rotated_coords[0,:,:]=input_coords.T


for i in range(1,n_side):
    rotated_coords[i,:,:]=rotate_side(n_side,input_coords)

    plt.plot(rotated_coords[i,0,:],rotated_coords[i,1,:],marker='x')
    i_side+=1
    input_coords=rotated_coords[i,:,:].T
plt.plot(rotated_coords[0,0,:],rotated_coords[0,1,:],marker='x')
plt.show()

rotated_coords_translated=rotated_coords+box_centre
for i in range(0,n_side):
    

    plt.plot(rotated_coords_translated[i,0,:],rotated_coords_translated[i,1,:],marker='x')

conv_r_t_coords=np.transpose(rotated_coords_translated,(0,2,1))
print(conv_r_t_coords)

#%% translate shape into chain
# this works for equilateral triangle 
number_of_units=3
coord_chain=np.zeros((number_of_units,n_side,2,3))
#s = 2r Sin(180/n)
side_length=np.abs(2*rad_circle*np.sin(np.pi/n_side))
for i in range(number_of_units):
    xshift=side_length*(i+1)

    coord_chain[i,:,:,:]=rotated_coords_translated+np.array([[xshift,xshift,xshift],
                                                             [0,0,0]])
    



for i in range(number_of_units):
    for j in range(n_side):


        plt.plot(coord_chain[i,j,0,:],coord_chain[i,j,1,:])
plt.axis('equal')
plt.plot()




#%% connected at a centre point 

number_of_units=3 
coord_chain=np.zeros((number_of_units,n_side,2,3))
#s = 2r Sin(180/n)
side_length=np.abs(2*rad_circle*np.sin(np.pi/n_side))
xshift=side_length
coord_chain[0,:,:,:]=rotated_coords_translated
xshift=side_length


                                                         
coord_chain[1,:,:,:]=rotated_coords_translated+np.array([[xshift,xshift,xshift],
                                                             [0,0,0]])                                                   
xshift=side_length/2
yshift=-np.sqrt((side_length**2)-((0.5*side_length)**2))
    

coord_chain[2,:,:,:]=rotated_coords_translated+np.array([[xshift,xshift,xshift],
                                                             [yshift,yshift,yshift]])   



for i in range(3):
    for j in range(n_side):


        plt.plot(coord_chain[i,j,0,:],coord_chain[i,j,1,:])
plt.axis('equal')
plt.plot()



#%% check shape 



#%% create molecule file 

