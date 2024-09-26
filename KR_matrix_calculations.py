#%%
import numpy as np

M=np.array([[2,-1,0],
        [-1,1,0],
        [0,0,1]])
determinant= np.linalg.det(M)

Eigennumbers=np.linalg.eig(M)
# %%
V=Eigennumbers[1] # eigenvectors
Lambda=Eigennumbers[0] # eigenvalues
Lambda_matrix=np.array([[2.61803399,0,0],
        [0,0.38196601,0],
        [0,0,1]])

# test
LHS=np.matmul(M,V)
RHS=np.matmul(V,Lambda_matrix)
print(LHS)
print(RHS)

V_inv=np.linalg.inv(V)
print(V_inv)
# %%
