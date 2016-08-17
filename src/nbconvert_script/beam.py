
# coding: utf-8

# In[16]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/4.4.png" alt="Figure 4.4" align=left>

# In[17]:

# Global
var('E')
# Element
var('A_1 A_2')
var('I_1 I_2')
var('l_1 l_2')
var('w_1 w_2')
# Node
var('U_11 U_12 U_21 U_22 U_31 U_32')
var('R_11 R_12 R_21 R_22 R_31 R_32');


# In[18]:

# Predefined Value
E_val = 200e9
w_1_val = w_2_val = -25000
I_1_val = I_2_val = 1.186e-4
l_1_val, l_2_val = 5, 2.5


# In[19]:

# Global Vector & Matrix
K_1 = MatMul(E*I_1/l_1**3, Matrix([[12, 6*l_1, -12, 6*l_1], 
                                 [6*l_1, 4*l_1**2, -6*l_1, 2*l_1**2],
                                 [-12, -6*l_1, 12, -6*l_1],
                                 [6*l_1, 2*l_1**2, -6*l_1, 4*l_1**2]]))
K_2 = K_1.subs([(l_1, l_2), (I_1, I_2)])
K_G = zeros(6, 6)
K_G[0:4, 0:4] += K_1
K_G[2:6, 2:6] += K_2

F_1 = Matrix(4, 1, [w_1*l_1/2, w_1*l_1**2/12, w_1*l_1/2, -w_1*l_1**2/12])
F_2 = Matrix(4, 1, [w_2*l_2/2, w_2*l_2**2/12, w_2*l_2/2, -w_2*l_2**2/12])
F_G = zeros(6, 1)
F_G[0:4, 0] += F_1
F_G[2:6, 0] += F_2

U_1 = Matrix(4, 1, [U_11, U_12, U_21, U_22])
U_2 = Matrix(4, 1, [U_21, U_22, U_31, U_32])
U_G = Matrix(6, 1, [U_11, U_12, U_21, U_22, U_31, U_32])

R_1 = Matrix(4, 1, [R_11, R_12, R_21, R_22])
R_2 = Matrix(4, 1, [R_21, R_22, R_31, R_32])
R_G = Matrix(6, 1, [R_11, R_12, R_21, R_22, R_31, R_32])

print('>>> Global Stiffness Matrix:')
eqs_disp(['K_1', K_1], ['K_2', K_2])
eq_disp('[K_G]', K_G)
print('>>> Global External Force Vector:')
eqs_disp(['F_1', F_1], ['F_2', F_2])
eq_disp('[F_G]', F_G)
print('>>> Global Displacement Vector:')
eqs_disp(['U_1', U_1], ['U_2', U_2])
eq_disp('[U_G]', U_G)
print('>>> Global Reaction Force Vector:')
eqs_disp(['R_1', R_1], ['R_2', R_2])
eq_disp('[R_G]', R_G)


# In[20]:

# Global Stiffness Matrix (Reaction Force Eliminated)
K_G_val = K_G.subs([(I_1, I_1_val), (I_2, I_2_val),
                    (E, E_val),
                    (l_1, l_1_val), (l_2, l_2_val)])

K_GM_val = Matrix(K_G_val)
K_GM_val[0:3, 0:3] = eye(3)
K_GM_val[0:3, 3:6] = zeros(3, 3)
K_GM_val[3:6, 0:3] = zeros(3, 3)

print('>>> Global Stiffness Matrix (Reaction Force Eliminated):')
eq_disp('[K_{GM}]', K_GM_val)


# In[21]:

# Solve
F_G_val = F_G.subs([(w_1, w_1_val), (w_2, w_2_val), (l_1, l_1_val), (l_2, l_2_val)])

F_GM_val = Matrix(F_G_val)
F_GM_val[0:3, 0] = Matrix(3, 1, [0, 0, 0])

U_G_val = K_GM_val ** -1 * F_GM_val
R_G_val = K_G_val * U_G_val - F_G_val

print('>>> Global External Force Vector:')
eq_disp('\{F_G\}', F_G_val.evalf())
eq_disp('\{F_{GM}\}', F_GM_val.evalf())
print('>>> Global Displacement Vector:')
eq_disp('\{U_G\}', U_G_val)
print('>>> Global Reaction Force Vector:')
eq_disp('\{R_G\}', R_G_val)

