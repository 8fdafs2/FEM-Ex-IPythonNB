
# coding: utf-8

# In[6]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/1.1.png" alt="Figure 1.1" align=left>
# <img src="image/1.2.png" alt="Figure 1.2" align=left>
# <img src="image/1.4.png" alt="Figure 1.4" align=left>
# <img src="image/1.5.png" alt="Figure 1.5" align=left>

# In[7]:

# Global
var('w_1 w_2 L P t E')
# Element
var('l_1 l_2 l_3 l_4')
var('k_1 k_2 k_3 k_4')
var('epsilon_1 epsilon_2 epsilon_3 epsilon_4')
# Node
var('A_1 A_2 A_3 A_4 A_5') 
var('u_1 u_2 u_3 u_4 u_5')
var('f_1 f_2 f_3 f_4 f_5')
var('R_1 R_2 R_3 R_4 R_5')
var('P_1 P_2 P_3 P_4 P_5');


# In[8]:

# Predefined Value
E_val = 10.4e6
w_1_val = 2
w_2_val = 1
t_val = 0.125
L_val = 10
P_val = 1e3


# In[9]:

K_1 = Matrix([[k_1, -k_1], [-k_1, k_1]]); 
K_2 = Matrix([[k_2, -k_2], [-k_2, k_2]]); 
K_3 = Matrix([[k_3, -k_3], [-k_3, k_3]]); 
K_4 = Matrix([[k_4, -k_4], [-k_4, k_4]]); 

# Direct Method

print('>>> Element Stiffness Matrix by Direct Method:')
eqs_disp( [ Matrix([f_1, f_2]), MatMul(K_1, Matrix([u_1, u_2])) ],
          [ Matrix([f_2, f_3]), MatMul(K_2, Matrix([u_2, u_3])) ],
          [ Matrix([f_3, f_4]), MatMul(K_3, Matrix([u_3, u_4])) ],
          [ Matrix([f_4, f_5]), MatMul(K_4, Matrix([u_4, u_5])) ],
         inline = False)

# Energy Method

Lambda_1 = Function('Lambda_1')(u_1, u_2)
Lambda_2 = Function('Lambda_2')(u_2, u_3)
Lambda_3 = Function('Lambda_3')(u_3, u_4)
Lambda_4 = Function('Lambda_4')(u_4, u_5)

print('>>> Element Stiffness Matrix by Energy Method:')
eqs_disp( [ Matrix([Derivative(Lambda_1, u_1), Derivative(Lambda_1, u_2)]), MatMul(K_1, Matrix([u_1, u_2])) ],
          [ Matrix([Derivative(Lambda_2, u_2), Derivative(Lambda_2, u_3)]), MatMul(K_2, Matrix([u_2, u_3])) ],
          [ Matrix([Derivative(Lambda_3, u_3), Derivative(Lambda_3, u_4)]), MatMul(K_3, Matrix([u_3, u_4])) ],
          [ Matrix([Derivative(Lambda_4, u_4), Derivative(Lambda_4, u_5)]), MatMul(K_4, Matrix([u_4, u_5])) ],
          inline = False)

K = zeros(5,5)
K[0:2, 0:2] += K_1
K[1:3, 1:3] += K_2
K[2:4, 2:4] += K_3
K[3:5, 3:5] += K_4

u = Matrix([u_1, u_2, u_3, u_4, u_5])

R_1, R_2, R_3, R_4, R_5 = R_1, 0, 0, 0, 0
R = Matrix([R_1, R_2, R_3, R_4, R_5])

P_1, P_2, P_3, P_4, P_5 = 0, 0, 0, 0, P_5
F = Matrix([P_1, P_2, P_3, P_4, P_5])

K_m = Matrix(K)
K_m[0, 0] = eye(1)
K_m[1:5, 0] = zeros(4, 1)
K_m[0, 1:5] = zeros(1, 4)

print('>>> Equilibrium Equation:')
eq_disp('\{R\}', '[K]\{u\} - \{F\}')
eq_disp(R, MatMul(K, u), -F)

print('>>> Equilibrium Equation (Reaction Force eliminated):')
eq_disp('[K_m]\{u\}', '\{F\}')
eq_disp(MatMul(K_m, u), F)


# In[10]:

k_1 = (A_2 + A_1) * E / (2 * l_1)
k_2 = (A_3 + A_2) * E / (2 * l_2)
k_3 = (A_4 + A_3) * E / (2 * l_3)
k_4 = (A_5 + A_4) * E / (2 * l_4)

l_1_val = l_2_val = l_3_val = l_4_val = L_val / 4

A_1_val = t_val * w_1_val
A_2_val = t_val * (w_1_val - (w_1_val - w_2_val) / 4 )
A_3_val = t_val * (w_1_val - (w_1_val - w_2_val) / 4 * 2 )
A_4_val = t_val * (w_1_val - (w_1_val - w_2_val) / 4 * 3 )
A_5_val = t_val * w_2_val

k_1_val = k_1.subs([(E, E_val), (A_1, A_1_val), (A_2, A_2_val), (l_1, l_1_val)])
k_2_val = k_2.subs([(E, E_val), (A_2, A_2_val), (A_3, A_3_val), (l_2, l_2_val)])
k_3_val = k_3.subs([(E, E_val), (A_3, A_3_val), (A_4, A_4_val), (l_3, l_3_val)])
k_4_val = k_4.subs([(E, E_val), (A_4, A_4_val), (A_5, A_5_val), (l_4, l_4_val)])

print('>>> Element Stiffness:')
eqs_disp(('k_1', k_1), ('k_2', k_2), ('k_3', k_3), ('k_4', k_4))
eqs_disp(('k_1', k_1_val), ('k_2', k_2_val), ('k_3', k_3_val), ('k_4', k_4_val))


# In[11]:

K_m_val = K_m.subs([('k_1', k_1_val), ('k_2', k_2_val), ('k_3', k_3_val), ('k_4', k_4_val)])

F_val = F.subs(P_5, P_val)

u_val = K_m_val**-1 * F_val

K_val = K.subs([('k_1', k_1_val), ('k_2', k_2_val), ('k_3', k_3_val), ('k_4', k_4_val)])

R_val = K_val * u_val - F_val

print('>>> Displacement Vector:')
eq_disp('\{u\}', 
        '[K_m]^{-1}\{F\}', True, 
        MatPow(K_m_val, -1) * F_val, True,
        u_val)

print('>>> Reaction Force Vector:')
eq_disp('\{R\}', 
        '[K]\{u\} - \{F\}', True,
        MatMul(K_val, u_val) - F_val, True,
        R_val)


# In[12]:

sigma_1 = E * epsilon_1
sigma_2 = E * epsilon_2
sigma_3 = E * epsilon_3
sigma_4 = E * epsilon_4

u_1_val, u_2_val, u_3_val, u_4_val, u_5_val = u_val

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


# In[ ]:



