
# coding: utf-8

# In[30]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/4.7.png" alt="Figure 4.7" align=left>

# In[31]:

# Global
var('E')
var('l_e')
var('u_i u_j sigma_i sigma_j epsilon_i epsilon_j')
# Element
var('A_1 A_2 A_3 A_4')
var('l_1 l_2 l_3 l_4')
# Node
var('F_1 F_2 F_3 F_4 F_5')
var('U_1 U_2 U_3 U_4 U_5')
var('R_1 R_2 R_3 R_4 R_5')
var('epsilon_1 epsilon_2 epsilon_3 epsilon_4');


# In[32]:

# Predefined Value
E_val = 29e6
A_1_val = A_2_val = A_3_val = A_4_val = 39.7
l_1_val = l_2_val = l_3_val = l_4_val = 15 * 12
F_1_val = 0; F_2_val = F_3_val = F_4_val = 25000 * 2; F_5_val = 30000 * 2


# In[33]:

# Global Stiffness Matrix
K_1 = MatMul(A_1 * E / l_1, Matrix([[1, -1],[-1, 1]]))
K_2 = MatMul(A_2 * E / l_2, Matrix([[1, -1],[-1, 1]]))
K_3 = MatMul(A_3 * E / l_3, Matrix([[1, -1],[-1, 1]]))
K_4 = MatMul(A_4 * E / l_4, Matrix([[1, -1],[-1, 1]]))

K_G = zeros(5, 5)
K_G[0:2, 0:2] += K_1
K_G[1:3, 1:3] += K_2
K_G[2:4, 2:4] += K_3
K_G[3:5, 3:5] += K_4

print('>>> Global Stiffness Matrix:')
eqs_disp(['[K_1]', K_1], ['[K_2]', K_2], ['[K_3]', K_3], ['[K_4]', K_4])
eq_disp('[K_G]', K_G)


# In[34]:

# Global Stiffness Matrix (Reaction Force Eliminated)
R = Matrix([R_1, R_2, R_3, R_4, R_5])
U = Matrix([U_1, U_2, U_3, U_4, U_5])
F = Matrix([F_1, F_2, F_3, F_4, F_5])

K_GM = Matrix(K_G)
K_GM[0, 0:5] = Matrix(1, 5, [1,0,0,0,0])

print('>>> Global Stiffness Matrix (Reaction Force Eliminated):')
eq_disp('[R]', '[K_G]\{U\}-\{F\}')
eq_disp(R, MatMul(K_G, U) - F)
print('~~~~~~~~~~~')
eq_disp('[K_{GM}]\{U\}', '\{F\}')
eq_disp(MatMul(K_GM, U), F)


# In[35]:

# Global Stiffness Matrix (Reaction Force Eliminated)
K_G_val = K_G.subs([(A_1, A_1_val), (A_2, A_2_val), (A_3, A_3_val), (A_4, A_4_val),
                    (E, E_val),
                    (l_1, l_1_val), (l_2, l_2_val), (l_3, l_3_val), (l_4, l_4_val),
                   ])

K_GM_val = Matrix(K_G_val)

K_GM_val[0, 0:5] = Matrix(1, 5, [1,0,0,0,0])
K_GM_val[1:5, 0] = zeros(4, 1)

print('>>> Global Stiffness Matrix (Reaction Force Eliminated):')
eq_disp('[K_{GM}]', K_GM_val)


# In[36]:

# Solve Displacement
F_val = F.subs([(F_1, F_1_val), (F_2, F_2_val), (F_3, F_3_val), (F_4, F_4_val), (F_5, F_5_val)])
U_val = K_GM_val ** -1 * F_val
R_val = K_G_val * U_val - F_val

print('>>> Global External Force Vector:')
eq_disp('\{F\}', F_val)
print('>>> Global Displacement Vector:')
eq_disp('\{U\}', U_val)
print('>>> Global Reaction Force Vector:')
eq_disp('\{R\}', R_val)


# In[37]:

# Solve Strain & Stress
sigma_1 = E * epsilon_1
sigma_2 = E * epsilon_2
sigma_3 = E * epsilon_3
sigma_4 = E * epsilon_4

u_1_val, u_2_val, u_3_val, u_4_val, u_5_val = U_val

epsilon_1 = (u_2 - u_1)/l_1
epsilon_2 = (u_3 - u_2)/l_2
epsilon_3 = (u_4 - u_3)/l_3
epsilon_4 = (u_5 - u_4)/l_4

epsilon_1_val = epsilon_1.subs([(u_1, u_1_val), (u_2, u_2_val), (l_1, l_1_val)])
epsilon_2_val = epsilon_2.subs([(u_2, u_2_val), (u_3, u_3_val), (l_2, l_2_val)])
epsilon_3_val = epsilon_3.subs([(u_3, u_3_val), (u_4, u_4_val), (l_3, l_3_val)])
epsilon_4_val = epsilon_4.subs([(u_4, u_4_val), (u_5, u_5_val), (l_4, l_4_val)])

epsilon = Matrix([epsilon_1, epsilon_2, epsilon_3, epsilon_4])
epsilon_val = Matrix([epsilon_1_val, epsilon_2_val, epsilon_3_val, epsilon_4_val])

sigma_1_val = sigma_1.subs([(E, E_val), ('epsilon_1', epsilon_1_val)])
sigma_2_val = sigma_2.subs([(E, E_val), ('epsilon_2', epsilon_2_val)])
sigma_3_val = sigma_3.subs([(E, E_val), ('epsilon_3', epsilon_3_val)])
sigma_4_val = sigma_4.subs([(E, E_val), ('epsilon_4', epsilon_4_val)])

sigma = Matrix([sigma_1, sigma_2, sigma_3, sigma_4])
sigma_val = Matrix([sigma_1_val, sigma_2_val, sigma_3_val, sigma_4_val])

print('>>> Strain Vector:')
eq_disp('\{\epsilon\}', epsilon, epsilon_val)
print('>>> Stress Vector:')
eq_disp('\{\sigma\}', sigma, sigma_val)


# In[38]:

# Solve on Certain Point (Not Node)
pos = 33
e = 3
i = 3
j = 4
y_val = 3
l_e_val = l_3_val
u_i_val = u_3_val
u_j_val = u_4_val
sigma_i_val = sigma_3_val
sigma_j_val = sigma_4_val
epsilon_i_val = epsilon_3_val
epsilon_j_val = epsilon_4_val

u_e = 1 / l_e * (u_i * (l_e - y) + u_j * y)
sigma_e = 1 / l_e * (sigma_i * (l_e - y) + sigma_j * y)
epsilon_e = 1 / l_e * (epsilon_i * (l_e - y) + epsilon_j * y)

u_e_val = u_e.subs([(l_e, l_e_val), (y, y_val), 
                    (u_i, u_i_val), (u_j, u_j_val)])

sigma_e_val = sigma_e.subs([(l_e, l_e_val), (y, y_val), 
                    (sigma_i, sigma_i_val), (sigma_j, sigma_j_val)])

epsilon_e_val = epsilon_e.subs([(l_e, l_e_val), (y, y_val), 
                    (epsilon_i, epsilon_i_val), (epsilon_j, epsilon_j_val)])

print('>>> Solve on Point of Y=33:')
eq_disp('u_e(y)', r'\frac{1}{l_e}(u_i(l_e-y)+u_jy)', u_e_val)
eq_disp(r'\sigma_e(y)', r'\frac{1}{l_e}(\sigma_i(l_e-y)+\sigma_jy)', sigma_e_val)
eq_disp(r'\epsilon_e(y)', r'\frac{1}{l_e}(\epsilon_i(l_e-y)+\epsilon_jy)', epsilon_e_val)

