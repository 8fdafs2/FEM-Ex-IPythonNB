
# coding: utf-8

# In[1]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/1.1.png" alt="Figure 1.1" align=left>

# In[2]:

# Global
var('w_1 w_2 L P t E')
var('A, y, u, du_dy, epsilon')
# Element
var('l_1 l_2 l_3 l_4') 
var('epsilon_1 epsilon_2 epsilon_3 epsilon_4')
var('sigma_1 sigma_2 sigma_3 sigma_4')
# Node
var('u_1 u_2 u_3 u_4 u_5')


# In[5]:

# Predefined Value
E_val = 10.4e6
w_1_val = 2
w_2_val = 1
t_val = 0.125
L_val = 10
P_5_val = 1e3


# In[3]:

A = (w_1 + (w_2-w_1)/L*y)*t
du_dy = P / (E * A)

print('>>> Equilibrium Equation:')
eqs_disp([r'{\epsilon}', r'\frac{du(y)}{dy}'], [r'A(y)', A])
eq_disp(r'P-E{\epsilon}A(y)', r'P-EA(y)\frac{du(y)}{dy}', 0)
eq_disp(r'\frac{du(y)}{dy}', r'\frac{P}{EA(y)}', du_dy)


# In[4]:

u = Integral(du_dy, (y, 0, y))
u_ = simplify(u.doit())

print('>>> Diplacement: (obtained by integration)')

eq_disp('u(y)', u, True, 
        u_)

u = u_


# In[6]:

l_1_val = l_2_val = l_3_val = l_4_val = L_val / 4

u = u.subs([(E, E_val), (L, L_val), (P, P_5_val), (w_1, w_1_val), (w_2, w_2_val), (t, t_val)])

u_1_val = u.subs(y, 0)
u_2_val = u.subs(y, l_1_val)
u_3_val = u.subs(y, l_1_val + l_2_val)
u_4_val = u.subs(y, l_1_val + l_2_val + l_3_val)
u_5_val = u.subs(y, l_1_val + l_2_val + l_3_val + l_4_val)

u_val = Matrix([u_1_val, u_2_val, u_3_val, u_4_val, u_5_val])

print('>>> Diplacement Vector:')
eq_disp(r'\{u\}', u_val)


# In[7]:

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



