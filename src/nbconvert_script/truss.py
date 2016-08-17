
# coding: utf-8

# In[35]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <table style="border-style: hidden; border-width: 0px; font-size: 12px; border-collapse: collapse; border-color: gray; background-color: #F5F5F5">
#     <tr style="border-style: hidden;">
#         <td style="border-style: hidden;" colspan=2 align=center>
#             <img src="image/2.5.png" alt="Figure 2.5">
#             <center><b>Figure 2.5</b></center>
#         </td>
#     </tr>
#     <tr style="border-style: hidden;">
#         <td style="border-style: hidden;" align=center>
#             <img src="image/2.6.png" alt="Figure 2.6">
#             <center><b>Figure 2.6</b></center>
#         </td>
#         <td style="border-style: hidden;" align=center>
#             <img src="image/2.7.png" alt="Figure 2.7">
#             <center><b>Figure 2.7</b></center>
#         </td>
#     </tr>
#     <tr style="border-style: hidden;">
#         <td style="border-style: hidden;" align=center>
#             <img src="image/2.8.png" alt="Figure 2.8">
#             <center><b>Figure 2.8</b></center>
#         </td>
#         <td style="border-style: hidden;" align=center>
#             <img src="image/2.9.png" alt="Figure 2.9">
#             <center><b>Figure 2.9</b></center>
#         </td>
#     </tr>
#     <tr style="border-style: hidden;">
#         <td style="border-style: hidden;" colspan=2 align=center>
#             <img src="image/2.10.png" alt="Figure 2.10">
#             <center><b>Figure 2.10</b></center>
#         </td>
#     </tr>
# </table>

# In[36]:

# Element
var('A_1 A_2 A_3 A_4 A_5 A_6')
var('L_1 L_2 L_3 L_4 L_5 L_6')
var('E_1 E_2 E_3 E_4 E_5 E_6')
var('k_1 k_2 k_3 k_4 k_5 k_6')
var('theta_1 theta_2 theta_3 theta_4 theta_5 theta_6')
var('epsilon_1 epsilon_2 epsilon_3 epsilon_4 epsilon_5 epsilon_6')
# Node
var('U_1x U_1y U_2x U_2y U_3x U_3y U_4x U_4y U_5x U_5y')
var('u_1x u_1y u_2x u_2y u_3x u_3y u_4x u_4y u_5x u_5y')
var('F_1x F_1y F_2x F_2y F_3x F_3y F_4x F_4y F_5x F_5y')
var('f_1x f_1y f_2x f_2y f_3x f_3y f_4x f_4y f_5x f_5y')
var('R_1x R_1y R_2x R_2y R_3x R_3y R_4x R_4y R_5x R_5y');


# In[37]:

# Predefined Value
U_1x_val = U_1y_val = 0
U_3x_val = U_3y_val = 0
F_4y_val = F_5y_val = -500

theta_1_val, theta_2_val, theta_3_val, theta_4_val, theta_5_val, theta_6_val = 0, rad(135), 0, rad(90), rad(45), 0
E_1_val = E_2_val = E_3_val = E_4_val = E_5_val = E_6_val = 1.9e6
A_1_val = A_2_val = A_3_val = A_4_val = A_5_val = A_6_val = 8
L_1_val = L_3_val = L_4_val = L_6_val = 36
L_2_val = L_5_val = 50.9


# In[38]:

# Vectors in Local Coord. Sys. & Global Coord. Sys.
U_1 = Matrix([U_1x, U_1y, U_2x, U_2y]); u_1 = Matrix([u_1x, u_1y, u_2x, u_2y])
F_1 = Matrix([F_1x, F_1y, F_2x, F_2y]); f_1 = Matrix([f_1x, f_1y, f_2x, f_2y])
U_2 = Matrix([U_2x, U_2y, U_3x, U_3y]); u_2 = Matrix([u_2x, u_2y, u_3x, u_3y])
F_2 = Matrix([F_2x, F_2y, F_3x, F_3y]); f_2 = Matrix([f_2x, f_2y, f_3x, f_3y])
U_3 = Matrix([U_3x, U_3y, U_4x, U_4y]); u_3 = Matrix([u_3x, u_3y, u_4x, u_4y])
F_3 = Matrix([F_3x, F_3y, F_4x, F_4y]); f_3 = Matrix([f_3x, f_3y, f_4x, f_4y])
U_4 = Matrix([U_2x, U_2y, U_4x, U_4y]); u_4 = Matrix([u_2x, u_2y, u_4x, u_4y])
F_4 = Matrix([F_2x, F_2y, F_4x, F_4y]); f_4 = Matrix([f_2x, f_2y, f_4x, f_4y])
U_5 = Matrix([U_2x, U_2y, U_5x, U_5y]); u_5 = Matrix([u_2x, u_2y, u_5x, u_5y])
F_5 = Matrix([F_2x, F_2y, F_5x, F_5y]); f_5 = Matrix([f_2x, f_2y, f_5x, f_5y])
U_6 = Matrix([U_4x, U_4y, U_5x, U_5y]); u_6 = Matrix([u_4x, u_4y, u_5x, u_5y])
F_6 = Matrix([F_4x, F_4y, F_5x, F_5y]); f_6 = Matrix([f_4x, f_4y, f_5x, f_5y])


# In[39]:

# Transition Matrix for Converting Vector from Local Coord. Sys. to Global one
T_1 = Matrix([[cos(theta_1), -sin(theta_1), 0, 0], 
              [sin(theta_1),  cos(theta_1), 0, 0],
              [0, 0, cos(theta_1), -sin(theta_1)],
              [0, 0, sin(theta_1),  cos(theta_1)]])

T_2 = T_1.subs(theta_1, theta_2)
T_3 = T_1.subs(theta_1, theta_3)
T_4 = T_1.subs(theta_1, theta_4)
T_5 = T_1.subs(theta_1, theta_5)
T_6 = T_1.subs(theta_1, theta_6)

print('>>> Converting Displacement Vector from Local Coord. Sys. to Global one:')
eq_disp('\{U_1\}', '[T_1]\{u_1\}')
print('------------------')
eq_disp(U_1, MatMul(T_1, u_1))
print('\n'.join(list('...')))
eq_disp('\{U_6\}', '[T_6]\{u_6\}')
print('------------------')
eq_disp(U_6, MatMul(T_6, u_6))
print('\n')
print('>>> Converting Load Vector from Local Coord. Sys. to Global one:')
eq_disp('\{F_1\}', '[T_1]\{f_1\}')
print('------------------')
eq_disp(U_1, MatMul(T_1, u_1))
print('\n'.join(list('...')))
eq_disp('\{F_6\}', '[T_6]\{f_6\}')
print('------------------')
eq_disp(U_1, MatMul(T_1, u_1))


# In[40]:

# Stiffness Matrix in Local Coord. Sys.
K_1 = Matrix([[k_1,  0, -k_1, 0], 
              [0,    0, 0,    0],
              [-k_1, 0, k_1,  0],
              [0,    0, 0,    0]])

K_2 = K_1.subs(k_1, k_2)
K_3 = K_1.subs(k_1, k_3)
K_4 = K_1.subs(k_1, k_4)
K_5 = K_1.subs(k_1, k_5)
K_6 = K_1.subs(k_1, k_6)

print('>>> Stiffness Matrix in Local Coord. Sys.:')
eq_disp('\{f_1\}', '[K_1]\{u_1\}')
print('------------------')
eq_disp(f_1, MatMul(K_1, u_1))
print('\n'.join(list('...')))
eq_disp('\{f_6\}', '[K_6]\{u_6\}')
print('------------------')
eq_disp(f_6, MatMul(K_6, u_6))


# In[41]:

eq_disp('\{f\}', '[K]\{u\}')
eq_disp('[T]^{-1}\{F\}', '[K][T]^{-1}\{U\}')
eq_disp('\{F\}', '[T][K][T]^{-1}\{U\}')
eq_disp('\{F\}', '[K_G]\{U\}')


# In[42]:

# Stiffness Matrix in Global Coord. Sys.
K_G1 = simplify(T_1 * K_1 * T_1**-1)
K_G2 = K_G1.subs([(k_1, k_2), (theta_1, theta_2)])
K_G3 = K_G1.subs([(k_1, k_3), (theta_1, theta_3)])
K_G4 = K_G1.subs([(k_1, k_4), (theta_1, theta_4)])
K_G5 = K_G1.subs([(k_1, k_5), (theta_1, theta_5)])
K_G6 = K_G1.subs([(k_1, k_6), (theta_1, theta_6)])
# K_G2 = simplify(T_2 * K_2 * T_2**-1)
# K_G3 = simplify(T_3 * K_3 * T_3**-1)
# K_G4 = simplify(T_4 * K_4 * T_4**-1)
# K_G5 = simplify(T_5 * K_5 * T_5**-1)
# K_G6 = simplify(T_6 * K_6 * T_6**-1)

print('>>> Stiffness Matrix in Global Coord. Sys.:')
eq_disp('\{F_1\}', '[K_{G1}]\{U_1\}')
print('------------------')
eq_disp(F_1, MatMul(K_G1, U_1))
print('\n'.join(list('...')))
eq_disp('\{F_6\}', '[K_{G6}]\{U_6\}')
print('------------------')
eq_disp(F_6, MatMul(K_G6, U_6))


# In[43]:

# Global Stiffness Matrix
K_G = zeros(10,10)
K_G[0:4,  0:4] += K_G1
K_G[2:6,  2:6] += K_G2
K_G[4:8,  4:8] += K_G3
K_G[2:4,  2:4] += K_G4[0:2, 0:2]
K_G[2:4,  6:8] += K_G4[0:2, 2:4]
K_G[6:8,  2:4] += K_G4[2:4, 0:2]
K_G[6:8,  6:8] += K_G4[2:4, 2:4]
K_G[2:4,  2:4] += K_G5[0:2, 0:2]
K_G[2:4, 8:10] += K_G5[0:2, 2:4]
K_G[8:10, 2:4] += K_G5[2:4, 0:2]
K_G[8:10,8:10] += K_G5[2:4, 2:4]
K_G[6:10, 6:10] += K_G6

K_G_ = K_G.subs([(theta_1, theta_1_val), (theta_2, theta_2_val), (theta_3, theta_3_val), 
                 (theta_4, theta_4_val), (theta_5, theta_5_val), (theta_6, theta_6_val)])

print('>>> Global Stiffness Matrix:')
eq_disp('[K_G]', True,
        K_G, True,
        K_G_)

K_G = K_G_


# In[44]:

# Global Stiffness Matrix (Reaction Force eliminated)
K_GM = Matrix(K_G)

K_GM[0:2, 0:2] = eye(2)
K_GM[0:2, 2:10] = zeros(2, 8)
K_GM[2:10, 0:2] = zeros(8, 2)

K_GM[4:6, 4:6] = eye(2)
K_GM[0:4, 4:6] = zeros(4, 2)
K_GM[4:6, 0:4] = zeros(2, 4)
K_GM[6:10, 4:6] = zeros(4, 2)
K_GM[4:6, 6:10] = zeros(2, 4)

U = Matrix([U_1x, U_1y, U_2x, U_2y, U_3x, U_3y, U_4x, U_4y, U_5x, U_5y])
R = Matrix([R_1x, 0, 0, 0, R_3x, R_3y, 0, 0, 0, 0])
F = Matrix([0, 0, 0, 0, 0, 0, 0, F_4y, 0, F_5y])

print('>>> Equilibrium Equation:')
eq_disp('\{R\}', '[K_G]\{U\} - \{F\}')
print('------------------')
eq_disp(R, MatMul(K_G, U) - F)
print('\n')
print('>>> Global Stiffness Matrix (Reaction Force eliminated)')
eq_disp('[K_{GM}]\{U\}', '\{F\}')
print('------------------')
eq_disp(MatMul(K_GM, U), F)


# In[45]:

# Element Stiffness
k_1 = A_1 * E_1 / L_1
k_2 = A_2 * E_2 / L_2
k_3 = A_3 * E_3 / L_3
k_4 = A_4 * E_4 / L_4
k_5 = A_5 * E_5 / L_5
k_6 = A_6 * E_6 / L_6

k_1_val = k_1.subs([(A_1, A_1_val), (E_1, E_1_val), (L_1, L_1_val)])
k_2_val = k_2.subs([(A_2, A_2_val), (E_2, E_2_val), (L_2, L_2_val)])
k_3_val = k_3.subs([(A_3, A_3_val), (E_3, E_3_val), (L_3, L_3_val)])
k_4_val = k_4.subs([(A_4, A_4_val), (E_4, E_4_val), (L_4, L_4_val)])
k_5_val = k_5.subs([(A_5, A_5_val), (E_5, E_5_val), (L_5, L_5_val)])
k_6_val = k_6.subs([(A_6, A_6_val), (E_6, E_6_val), (L_6, L_6_val)])

print('>>> Element Stiffness:')
eqs_disp(('k_1', k_1), ('k_2', k_2), ('k_3', k_3), 
         ('k_4', k_4), ('k_5', k_5), ('k_6', k_6))

eqs_disp(('k_1', k_1_val), ('k_2', k_2_val), ('k_3', k_3_val))
eqs_disp(('k_4', k_4_val), ('k_5', k_5_val), ('k_6', k_6_val))


# In[46]:

# Displacement Solved in Regards to Global Coord. Sys.
K_GM_val = K_GM.subs([('k_1', k_1_val), ('k_2', k_2_val), ('k_3', k_3_val),
                      ('k_4', k_4_val), ('k_5', k_5_val), ('k_6', k_6_val)])

F_val = F.subs([('F_4y', F_4y_val), ('F_5y', F_5y_val)])

U_val = K_GM_val**-1 * F_val

K_G_val = K_G.subs([('k_1', k_1_val), ('k_2', k_2_val), ('k_3', k_3_val),
                    ('k_4', k_4_val), ('k_5', k_5_val), ('k_6', k_6_val)])

R_val = K_G_val * U_val - F_val

print('>>> Displacement Vector in Global Coord. Sys.:')
eq_disp('\{U\}', 
        '[K_{GM}]^{-1}\{F\}', True, 
        MatPow(K_GM_val, -1) * F_val)
eq_disp('\{U\}', U_val)
print('\n')
print('>>> Reaction Force Vector in Global Coord. Sys.:')
eq_disp('\{R\}', 
        '[K_G]\{U\} - \{F\}', True, 
        MatMul(K_G_val, U_val) - F_val)
eq_disp('\{R\}', R_val)


# In[49]:

# Displacement in Local Coord. Sys. Obtained by Using Transtion Matrix
print('>>> Displacement Vector in Local Coord. Sys.:')
eq_disp('\{u_1\}', '[T_1]^{-1}\{U_1\}')
print('------------------')
eq_disp(u_1, MatPow(T_1, -1) * U_1)
print('\n'.join(list('...')))
eq_disp('\{u_6\}', '[T_6]^{-1}\{U_6\}')
print('------------------')
eq_disp(u_6, MatPow(T_6, -1) * U_6)

U_1_val = Matrix(U_val[0:4])
T_1_val = T_1.subs(theta_1, theta_1_val)
u_1_val = T_1_val ** -1 * U_1_val

U_2_val = Matrix(U_val[2:6])
T_2_val = T_2.subs(theta_2, theta_2_val)
u_2_val = T_2_val ** -1 * U_2_val

U_3_val = Matrix(U_val[4:8])
T_3_val = T_3.subs(theta_3, theta_3_val)
u_3_val = T_3_val ** -1 * U_3_val

U_4_val = Matrix(U_val[2:4] + U_val[6:8])
T_4_val = T_4.subs(theta_4, theta_4_val)
u_4_val = T_4_val ** -1 * U_4_val

U_5_val = Matrix(U_val[2:4] + U_val[8:10])
T_5_val = T_5.subs(theta_5, theta_5_val)
u_5_val = T_5_val ** -1 * U_5_val

U_6_val = Matrix(U_val[6:10])
T_6_val = T_6.subs(theta_6, theta_6_val)
u_6_val = T_6_val ** -1 * U_6_val

print('\n')
eqs_disp(['\{u_1\}', u_1, u_1_val.evalf()], 
         ['\{u_2\}', u_2, u_2_val.evalf()],
         ['\{u_3\}', u_3, u_3_val.evalf()],
         ['\{u_4\}', u_4, u_4_val.evalf()],
         ['\{u_5\}', u_5, u_5_val.evalf()],
         ['\{u_6\}', u_6, u_6_val.evalf()])


# In[48]:

# Solve Strain & Stress on Basis of Displacement Solution
sigma_1 = E_1 * epsilon_1
sigma_2 = E_2 * epsilon_2
sigma_3 = E_3 * epsilon_3
sigma_4 = E_4 * epsilon_4
sigma_5 = E_5 * epsilon_5
sigma_6 = E_6 * epsilon_6

epsilon_1 = (u_2x - u_1x)/L_1
epsilon_2 = (u_3x - u_2x)/L_2
epsilon_3 = (u_4x - u_3x)/L_3
epsilon_4 = (u_4x - u_2x)/L_4
epsilon_5 = (u_5x - u_2x)/L_5
epsilon_6 = (u_5x - u_4x)/L_6

u_2x_val, u_1x_val = u_1_val[2], u_1_val[0]
epsilon_1_val = epsilon_1.subs([(u_2x, u_2x_val), (u_1x, u_1x_val), (L_1, L_1_val)])
u_3x_val, u_2x_val = u_2_val[2], u_2_val[0]
epsilon_2_val = epsilon_2.subs([(u_3x, u_3x_val), (u_2x, u_2x_val), (L_2, L_2_val)])
u_4x_val, u_3x_val = u_3_val[2], u_3_val[0]
epsilon_3_val = epsilon_3.subs([(u_4x, u_4x_val), (u_3x, u_3x_val), (L_3, L_3_val)])
u_4x_val, u_2x_val = u_4_val[2], u_4_val[0]
epsilon_4_val = epsilon_4.subs([(u_4x, u_4x_val), (u_2x, u_2x_val), (L_4, L_4_val)])
u_5x_val, u_2x_val = u_5_val[2], u_5_val[0]
epsilon_5_val = epsilon_5.subs([(u_5x, u_5x_val), (u_2x, u_2x_val), (L_5, L_5_val)])
u_5x_val, u_4x_val = u_6_val[2], u_6_val[0]
epsilon_6_val = epsilon_6.subs([(u_5x, u_5x_val), (u_4x, u_4x_val), (L_6, L_6_val)])

epsilon = Matrix([epsilon_1, epsilon_2, epsilon_3, epsilon_4, epsilon_5, epsilon_6])
epsilon_val = Matrix([epsilon_1_val, epsilon_2_val, epsilon_3_val, epsilon_4_val, epsilon_5_val, epsilon_6_val]).evalf()

print('>>> Strain Vector:')
eq_disp('\{\epsilon\}', epsilon, epsilon_val)

sigma_1_val = sigma_1.subs([(E_1, E_1_val), ('epsilon_1', epsilon_1_val)])
sigma_2_val = sigma_2.subs([(E_2, E_2_val), ('epsilon_2', epsilon_2_val)])
sigma_3_val = sigma_3.subs([(E_3, E_3_val), ('epsilon_3', epsilon_3_val)])
sigma_4_val = sigma_4.subs([(E_4, E_4_val), ('epsilon_4', epsilon_4_val)])
sigma_5_val = sigma_5.subs([(E_5, E_5_val), ('epsilon_5', epsilon_5_val)])
sigma_6_val = sigma_6.subs([(E_6, E_6_val), ('epsilon_6', epsilon_6_val)])

sigma = Matrix([sigma_1, sigma_2, sigma_3, sigma_4, sigma_5, sigma_6])
sigma_val = Matrix([sigma_1_val, sigma_2_val, sigma_3_val, sigma_4_val, sigma_5_val, sigma_6_val]).evalf()

print('>>> Stress Vector:')
eq_disp('\{\sigma\}', sigma, sigma_val)


# In[ ]:



