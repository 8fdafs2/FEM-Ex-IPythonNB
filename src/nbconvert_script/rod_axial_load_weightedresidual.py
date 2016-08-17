
# coding: utf-8

# In[2]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/1.1.png" alt="Figure 1.1" align=left>

# In[10]:

# Global
var('w_1 w_2 L P t E y')
var('c_1 c_2 c_3')
# Element
var('l_1 l_2 l_3 l_4') 
var('epsilon_1 epsilon_2 epsilon_3 epsilon_4')
# Node
var('u_1 u_2 u_3 u_4 u_5');


# In[4]:

# Predefined Value
E_val = 10.4e6
w_1_val = 2
w_2_val = 1
t_val = 0.125
L_val = 10
P_5_val = 1e3


# In[5]:

A = (w_1 + (w_2-w_1)/L*y)*t

print('>>> Equilibrium Equation:')
eqs_disp([r'{\epsilon}', r'\frac{du(y)}{dy}'], [r'A(y)', A])
eq_disp(r'P-E{\epsilon}A(y)', r'P-EA(y)\frac{du(y)}{dy}', 0)

u = c_1 * y + c_2 * y**2 + c_3 * y**3
du_dy = diff(u, y)

print('>>> Displacement assumed to respect the law of cubic polynomial:')
eq_disp(r'u(y)', u)
eq_disp(r'\frac{du}{dy}', du_dy)

Re = P - E * A * du_dy

print('>>> Residual:')
eq_disp(r'\Re', r'P-E{\epsilon}A(y)', Re)


# In[6]:

# Collocation Method

Eq_1_c = Re.subs(y, L/3)
Eq_2_c = Re.subs(y, L/3*2)
Eq_3_c = Re.subs(y, L)

Eq_1, Eq_2, Eq_3 = Eq_1_c, Eq_2_c, Eq_3_c

print('>>> Equation Set by Collocation Method: (In Use)')
eqs_disp(['Eq_1', Eq_1_c], ['Eq_2', Eq_2_c], ['Eq_3', Eq_3_c], inline=False)

# Subdomain Method

Eq_1_s = Integral(Re, [y, 0, L/3])
Eq_2_s = Integral(Re, [y, L/3, L/3*2])
Eq_3_s = Integral(Re, [y, L/3*2, L])

print('>>> Equation Set by Subdomain Method:')
eqs_disp(['Eq_1', Eq_1_s], ['Eq_2', Eq_2_s], ['Eq_3', Eq_3_s], inline=False)

# Galerkin Method

phi_1, phi_2, phi_3 = y, y**2, y**3

Eq_1_g = Integral(Re*phi_1, [y, 0, L])
Eq_2_g = Integral(Re*phi_2, [y, 0, L])
Eq_3_g = Integral(Re*phi_3, [y, 0, L])

print('>>> Equation Set by Galerkin Method:')
eqs_disp(['Eq_1', Eq_1_g], ['Eq_2', Eq_2_g], ['Eq_3', Eq_3_g], inline=False)

# LeastSquare Method

phi_1 = Derivative(Re, c_1)
phi_2 = Derivative(Re, c_2)
phi_3 = Derivative(Re, c_3)
Eq_1_l = Integral(Re*phi_1, [y, 0, L])
Eq_2_l = Integral(Re*phi_2, [y, 0, L])
Eq_3_l = Integral(Re*phi_3, [y, 0, L])

print('>>> Equation Set by LeastSquare Method:')
eqs_disp(['Eq_1', Eq_1_l], ['Eq_2', Eq_2_l], ['Eq_3', Eq_3_l], inline=False)


# In[7]:

Eq_1 = Eq_1.subs([(E, E_val), (P, P_5_val), (t, t_val), (L, L_val), (w_1, w_1_val), (w_2, w_2_val)]).doit()
Eq_2 = Eq_2.subs([(E, E_val), (P, P_5_val), (t, t_val), (L, L_val), (w_1, w_1_val), (w_2, w_2_val)]).doit()
Eq_3 = Eq_3.subs([(E, E_val), (P, P_5_val), (t, t_val), (L, L_val), (w_1, w_1_val), (w_2, w_2_val)]).doit()

sol = solve([Eq_1, Eq_2, Eq_3], [c_1, c_2, c_3])
c_1_val, c_2_val, c_3_val = sol[c_1], sol[c_2], sol[c_3]

print('>>> polynomial coefficient obtained by solving the Equation Set:')
eqs_disp([r'Eq_1', Eq_1, 0],
         [r'Eq_2', Eq_2, 0],
         [r'Eq_3', Eq_2, 0], inline = False)
eqs_disp([r'c_1', c_1_val],
         [r'c_2', c_2_val],
         [r'c_3', c_3_val], inline = False)


# In[8]:

l_1_val = l_2_val = l_3_val = l_4_val = L_val / 4

u = u.subs([(c_1, c_1_val), (c_2, c_2_val), (c_3, c_3_val)])

u_1_val = u.subs(y, 0)
u_2_val = u.subs(y, l_1_val)
u_3_val = u.subs(y, l_1_val + l_2_val)
u_4_val = u.subs(y, l_1_val + l_2_val + l_3_val)
u_5_val = u.subs(y, l_1_val + l_2_val + l_3_val + l_4_val)

u_val = Matrix([u_1_val, u_2_val, u_3_val, u_4_val, u_5_val])

print('>>> Diplacement Vector:')
eq_disp(r'\{u\}', u_val)


# In[9]:

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



