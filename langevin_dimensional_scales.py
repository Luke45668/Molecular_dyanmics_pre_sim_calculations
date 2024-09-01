#%%
import numpy as np

fluid_name="80%wt glycerol"
T_k=298 # fluid temp at density listed
real_fluid_viscosity=0.044 #Pa s 
K=np.array([40,50])
damp=0.035
k_b=1.38e-23
R_colloid=1e-7#1e-9
r_colloid_sim=0.25#0.25
mass_colloid_sim=5
sim_fluid_visc=mass_colloid_sim/(6*np.pi*r_colloid_sim*damp)
energy_scale= T_k*k_b
length_scale= R_colloid/r_colloid_sim
time_scale=real_fluid_viscosity* ((length_scale**3)/(sim_fluid_visc*energy_scale))
mass_scale=energy_scale*(time_scale**2)/(length_scale**2)

print("length scale",length_scale)
print("time scale", time_scale)
print("mass scale", mass_scale)
print("energy scale",energy_scale)

real_colloid_mass=5*mass_scale
erate=np.array([1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.175,0.15,0.125,0.1,0.08,
                0.06,0.04,0.02,0.01,0.005,0])
real_shear_rate=erate/time_scale
print("real shear rate",real_shear_rate)



# %%
