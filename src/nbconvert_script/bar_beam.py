
# coding: utf-8

# In[11]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/4.13.png" alt="Figure 4.13" align=left>
# <img src="image/4.14.png" alt="Figure 4.14" align=left>
# <img src="image/4.15.png" alt="Figure 4.15" align=left>

# In[12]:

# Global
var('E')
# Element
var('A_1 A_2')
var('I_1 I_2')
var('l_1 l_2')
var('theta_1 theta_2')
var('w_1 w_2')
# Node
var('U_11 U_12 U_13 U_21 U_22 U_23 U_31 U_32 U_33')
var('F_11 F_12 F_13 F_21 F_22 F_23 F_31 F_32 F_33')
var('R_11 R_12 R_13 R_21 R_22 R_23 R_31 R_32 R_33');


# In[13]:

# Predefined Value
E_val = 30e6
theta_1_val, theta_2_val = 0, rad(270)
A_1_val = A_2_val = 7.65
w_1_val, w_2_val = -800.0/12, 0
I_1_val = I_2_val = 204
l_1_val, l_2_val = 10*12, 9*12


# In[14]:

# Global Vector & Matrix
T_1 = eye(6)
T_1[0:2, 0:2] = Matrix([[cos(theta_1), sin(theta_1)], 
                        [-sin(theta_1), cos(theta_1)]])
T_1[3:5, 3:5] = Matrix([[cos(theta_1), sin(theta_1)], 
                        [-sin(theta_1), cos(theta_1)]])

T_2 = T_1.subs([(theta_1, theta_2)])
T_1_ = T_1
T_2_ = T_2
T_1 = T_1.subs([(theta_1, theta_1_val)]).evalf()
T_2 = T_2.subs([(theta_2, theta_2_val)]).evalf()

U_1 = Matrix(6, 1, [U_11, U_12, U_13, U_21, U_22, U_23])
U_2 = Matrix(6, 1, [U_21, U_22, U_23, U_31, U_32, U_33])
U_G = Matrix(9, 1, [U_11, U_12, U_13, U_21, U_22, U_23, U_31, U_32, U_33])

F_1 = Matrix(6, 1, [F_11, F_12, F_13, F_21, F_22, F_23])
F_2 = Matrix(6, 1, [F_21, F_22, F_23, F_31, F_32, F_33])
F_G = Matrix(9, 1, [F_11, F_12, F_13, F_21, F_22, F_23, F_31, F_32, F_33])
F_G_ = F_G
F_G = F_G.subs([(F_12, l_1*w_1/2), (F_13, l_1**2*w_1/12), 
                (F_21, 0), 
                (F_22, l_1*w_1/2), (F_23, -l_1**2*w_1/12)])

R_1 = Matrix(6, 1, [R_11, R_12, R_13, R_21, R_22, R_23])
R_2 = Matrix(6, 1, [R_21, R_22, R_23, R_31, R_32, R_33])
R_G = Matrix(9, 1, [R_11, R_12, R_13, R_21, R_22, R_23, R_31, R_32, R_33])


k_1 = Matrix([[A_1*E/l_1, 0, 0, -A_1*E/l_1, 0, 0], 
              [0, 12*E*I_1/l_1**3, 6*E*I_1/l_1**2, 0, -12*E*I_1/l_1**3, 6*E*I_1/l_1**2],
              [0, 6*E*I_1/l_1**2, 4*E*I_1/l_1, 0, -6*E*I_1/l_1**2, 2*E*I_1/l_1],
              [-A_1*E/l_1, 0, 0, A_1*E/l_1, 0, 0],
              [0, -12*E*I_1/l_1**3, -6*E*I_1/l_1**2, 0, 12*E*I_1/l_1**3, -6*E*I_1/l_1**2],
              [0, 6*E*I_1/l_1**2, 2*E*I_1/l_1, 0, -6*E*I_1/l_1**2, 4*E*I_1/l_1]])
k_2 = k_1.subs([(l_1, l_2), (I_1, I_2), (A_1, A_2)])
K_1 = T_1**-1 * k_1 * T_1
K_2 = T_2**-1 * k_2 * T_2
K_G = zeros(9, 9)
K_G[0:6, 0:6] += K_1
K_G[3:9, 3:9] += K_2

print('>>> Transition Matrix:')
eqs_disp(['T_1', T_1_, T_1], 
         ['T_2', T_2_, T_2])

print('>>> Global Stiffness Matrix:')
eqs_disp(['K_1', 'T_1^{-1}k_1T_1', MatMul(T_1**-1, k_1, T_1), K_1], 
         ['K_2', 'T_2^{-1}k_2T_2', MatMul(T_2**-1, k_2, T_2), K_2])
eq_disp('[K_G]', K_G)
print('>>> Global External Force Vector:')
eqs_disp(['F_1', F_1], 
         ['F_2', F_2])
eq_disp('[F_G]', F_G_, F_G)
print('>>> Global Displacement Vector:')
eqs_disp(['U_1', U_1], 
         ['U_2', U_2])
eq_disp('[U_G]', U_G)
print('>>> Global Reaction Force Vector:')
eqs_disp(['R_1', R_1], 
         ['R_2', R_2])
eq_disp('[R_G]', R_G)


# In[15]:

# Global Stiffness Matrix (Reaction Force Eliminated)
K_G_val = K_G.subs([(I_1, I_1_val), (I_2, I_2_val),
                    (E, E_val),
                    (l_1, l_1_val), (l_2, l_2_val),
                    (A_1, A_1_val), (A_2, A_2_val)])

K_GM_val = Matrix(K_G_val)
K_GM_val[0:3, 0:3] = eye(3)
K_GM_val[3:9, 0:3] = zeros(6, 3)
K_GM_val[0:3, 3:6] = zeros(3, 3)
K_GM_val[6:9, 3:6] = zeros(3, 3)
K_GM_val[0:6, 6:9] = zeros(6, 3)
K_GM_val[6:9, 6:9] = eye(3)

print('>>> Global Stiffness Matrix (Reaction Force Eliminated):')
eq_disp('[K_{GM}]', K_GM_val)


# In[16]:

# Solve
F_G_val = F_G.subs([(w_1, w_1_val),
                    (l_1, l_1_val)])

F_GM_val = Matrix(F_G_val)
F_GM_val[0:3, 0] = Matrix(3, 1, [0, 0, 0])
F_GM_val[6:9, 0] = Matrix(3, 1, [0, 0, 0])

U_G_val = K_GM_val ** -1 * F_GM_val
R_G_val = K_G_val * U_G_val - F_GM_val

print('>>> Global External Force Vector:')
eq_disp('\{F_G\}', F_G_val.evalf())
eq_disp('\{F_{GM}\}', F_GM_val.evalf())
print('>>> Global Displacement Vector:')
eq_disp('\{U_G\}', U_G_val)
print('>>> Global Reaction Force Vector:')
eq_disp('\{R_G\}', R_G_val)


# In[ ]:



