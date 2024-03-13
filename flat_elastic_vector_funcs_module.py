
import numpy as np

def bead_general_coordinates_eq_triangle(
    x,
    y,
    eqsl
    ):
    r_1= np.array([
        x,
        y
        ])
    r_2= np.array([
        x+eqsl,
        y
        ])
    r_3= np.array([
        x+(eqsl/2),
        y+(np.sqrt(eqsl)/2)*eqsl
        ])
    r_1_p=np.array([
        x+(eqsl/2),
        y
        ])
    r_2_p=np.array([
        x+(3*eqsl/4),
        y+(np.sqrt(eqsl)/4)*eqsl
        ])
    r_3_p=np.array([
         x+(eqsl/4),
         y+(np.sqrt(eqsl)/4)*eqsl
        ])
    coordinate_tuple=(
         r_1,
         r_2,
         r_3,
         r_1_p,
         r_2_p,
         r_3_p
        )

    return coordinate_tuple

def equilateral_coords_from_basis(first_coord, basis_vector_1,basis_vector_2,eqsl):
    r_1=first_coord
    
    r_2=(first_coord+ 
         basis_vector_1*eqsl
        )

    r_3=(first_coord+ 
        basis_vector_1*(eqsl/2) +
        basis_vector_2*(np.sqrt(eqsl)/2)*eqsl
        )
    r_1_p=(first_coord+
        basis_vector_1*(eqsl/2)
        )
    r_2_p=(first_coord +
        basis_vector_1*(3*eqsl/4)+ 
        basis_vector_2*(np.sqrt(eqsl)/4)*eqsl
        )
    r_3_p=(first_coord+
         basis_vector_1*(eqsl/4)+
         basis_vector_2*(np.sqrt(eqsl)/4)*eqsl
         )
    coordinate_tuple=(
         r_1,
         r_2,
         r_3,
         r_1_p,
         r_2_p,
         r_3_p
        )

    return coordinate_tuple










def point_on_plane(
          box_size_bar,
          box_size_index,
                   ):
      central_point=box_size_bar[box_size_index]*0.5
      p_0= np.array([
                   central_point,
                   central_point,
                   central_point
                   ])
      return p_0

def compute_vector_magnitude(vector):

    vector_mag= np.sqrt(vector[0]**2 +vector[1]**2 +vector[2]**2)

    return vector_mag

def compute_unit_vector(vector):

    vector_mag= np.sqrt(vector[0]**2 +vector[1]**2 +vector[2]**2)
    unit_vector=vector/vector_mag

    return unit_vector
    

def solve_plane_equation_for_z(
    normal,
    p_0,
    i, 
    j,
    coordinate_tuple_2d,
        ):
    
    z=-((normal[0]*(coordinate_tuple_2d[i][0]-p_0[0]) 
                 + normal[1]*(coordinate_tuple_2d[i][1]-p_0[1]))/normal[2])  + p_0[2]
    
    return z


def create_coordinates_string(coordinate_tuple_3d):
    coordinate_string=f"""1 {str(coordinate_tuple_3d[0][0])} {str(coordinate_tuple_3d[0][1])} {str(coordinate_tuple_3d[0][2])} \n2 {str(coordinate_tuple_3d[1][0])} {str(coordinate_tuple_3d[1][1])} {str(coordinate_tuple_3d[1][2])} \n3 {str(coordinate_tuple_3d[2][0])} {str(coordinate_tuple_3d[2][1])} {str(coordinate_tuple_3d[2][2])} \n4 {str(coordinate_tuple_3d[3][0])} {str(coordinate_tuple_3d[3][1])} {str(coordinate_tuple_3d[3][2])} \n5 {str(coordinate_tuple_3d[4][0])} {str(coordinate_tuple_3d[4][1])} {str(coordinate_tuple_3d[4][2])} \n6 {str(coordinate_tuple_3d[5][0])} {str(coordinate_tuple_3d[5][1])} {str(coordinate_tuple_3d[5][2])} """
    
    return coordinate_string


         

# compute unit rotation vector 


# from mol_file_overwriter import *

# #j=50
# for j in range(1):
#     coordinate_tuple_3d=()
#     box_size_index=4
#     p_0=point_on_plane(
#             box_size_bar,
#             box_size_index,
#                     )
#     coordinate_tuple_2d=bead_general_coordinates_eq_triangle(
#                                                             p_0[0],
#                                                             p_0[1],
#                                                             equilibrium_triangle_side_length
#                                                                 )



#     z=solve_plane_equation_for_z(
#             Rotation_vector[:,j],
#             p_0,
#             0,
#             j, 
#             coordinate_tuple_2d,
#                 )
#     first_coord = np.array([coordinate_tuple_2d[i][0],coordinate_tuple_2d[i][1],z])
#     basis_vector_1=(first_coord-p_0)/compute_vector_magnitude(first_coord-p_0)   
#     basis_vector_2= np.cross(Rotation_vector[:,j],(first_coord-p_0))/compute_vector_magnitude(np.cross(Rotation_vector[:,j],(first_coord-p_0)))
#     print(basis_vector_1)
#     print(basis_vector_2)


#     # now compute new coordinates in terms of basis 
    

#     coordinate_tuple_3d=equilateral_coords_from_basis(first_coord, basis_vector_1,basis_vector_2,equilibrium_triangle_side_length)



#     ell_1=coordinate_tuple_3d[1]-coordinate_tuple_3d[0]
#     ell_1_mag=compute_vector_magnitude(ell_1)
#     print(ell_1_mag)
#     ell_2=coordinate_tuple_3d[2]-coordinate_tuple_3d[0]
#     ell_2_mag=compute_vector_magnitude(ell_2)

    
#     print(compute_vector_magnitude(np.cross(ell_1,ell_2)))


#     fig = plt.figure()
#     ax = plt.axes(projection ='3d')

#     for i in range(len(coordinate_tuple_2d)):
#         x=coordinate_tuple_3d[i][0]
#         y=coordinate_tuple_3d[i][1]
#         z=coordinate_tuple_3d[i][2]
#         ax.scatter(x,y,z)
#     plt.show()

#     coord_string=create_coordinates_string(coordinate_tuple_3d)
#     name_particle='flatelastic_'+str(equilibrium_triangle_side_length)
#     rounded_normal= [ sgf.round(Rotation_vector[x,0],sigfigs=2) for x in range(3) ]


#     path_2_mol_template="/Volumes/Backup Plus 1/PhD_/Rouse Model simulations"\
#     "/Using LAMMPS imac/generic_molecule_files/flat_elastic_plate/"

#     Path_2_generic=path_2_mol_template
#     filename_generic="mol.elastic_plate_generic"
#     simulation_batch_folder="/Volumes/Backup Plus 1/PhD_/Rouse Model simulations"\
#     "/Using LAMMPS imac/generic_molecule_files/flat_elastic_plate/test_folder"

#     mol_file_coordinate_overwriter(Path_2_generic,
#                                    filename_generic,
#                                    simulation_batch_folder,
#                                    coord_string,
#                                    name_particle,
#                                    rounded_normal)



#     # now need to find template mol file 
#     # replace 
