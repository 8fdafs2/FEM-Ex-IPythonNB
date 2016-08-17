
# coding: utf-8

# In[13]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/4.12.png" alt="Figure 4.12" align=left>

# In[14]:

# Global
var('A E I l w theta')
var('F_i1 F_j1 F_i2 F_j2 F_i3 F_j3')
var('u_i1 u_j1 u_i2 u_j2 u_i3 u_j3');


# In[15]:

# Axial Stiffness Matrix
u_axial = Matrix([u_i1, u_j1])
F_axial = Matrix([F_i1, F_j1])
K_axial = MatMul(A * E / l, Matrix([[1, -1], [-1, 1]]))

print('>>> Axial Stffness Matrix:')
eq_disp(r'K_{axial}^{(e)}u_{axial}^{(e)}', 'F_{axial}^{(e)}')
print('~~~~~~~~~~~')
eq_disp('K_{axial}', K_axial)
eq_disp('F_{axial}', F_axial)
eq_disp('u_{axial}', u_axial)


# In[16]:

# XY Stiffness Matrix
u_xy = Matrix([u_i2, u_i3, u_j2, u_j3])
F_xy = Matrix([F_i2, F_i3, F_j2, F_j3])
F_xy_ = F_xy.subs([(F_i2, l*w/2), 
                   (F_i3, l**2*w/12), 
                   (F_j2, l*w/2), 
                   (F_j3, -l**2*w/12)])
K_xy = MatMul(E * I / l**3, Matrix([[12, 6*l, -12, 6*l], 
                                    [6*l, 4*l**2, -6*l, 2*l**2],
                                    [-12, -6*l, 12, -6*l],
                                    [6*l, 2*l**2, -6*l, 4*l**2]]))

print('>>> XY Stffness Matrix:')
eq_disp(r'K_{xy}^{(e)}u_{xy}^{(e)}', 'F_{xy}^{(e)}')
print('~~~~~~~~~~~')
eq_disp('K_{xy}', K_xy)
eq_disp('F_{xy}', F_xy, F_xy_)
eq_disp('u_{xy}', u_xy)


# In[17]:

# United Stiffness Matrix
u_united = Matrix([u_i1, u_i2, u_i3, u_j1, u_j2, u_j3])
F_united = Matrix([F_i1, F_i2, F_i3, F_j1, F_j2, F_j3])
F_united_ = F_united.subs([(F_i2, l*w/2), 
                           (F_i3, l**2*w/12), 
                           (F_j2, l*w/2), 
                           (F_j3, -l**2*w/12)])
K_united = zeros(6, 6)
K_united[0, 0] = K_axial[0, 0]
K_united[3, 0] = K_axial[1, 0]
K_united[0, 3] = K_axial[0, 1]
K_united[3, 3] = K_axial[1, 1]
K_united[1:3, 1:3] = K_xy[0:2, 0:2]
K_united[4:6, 1:3] = K_xy[2:4, 0:2]
K_united[1:3, 4:6] = K_xy[0:2, 2:4]
K_united[4:6, 4:6] = K_xy[2:4, 2:4]

print('>>> United Stffness Matrix:')
eq_disp('K_{united}', K_united)
eq_disp('F_{united}', F_united, F_united_)
eq_disp('u_{united}', u_united)


# In[18]:

# Transition Matrix
T_united = zeros(6, 6)
T_united[0:2, 0:2] = Matrix([[cos(theta), sin(theta)], 
                             [-sin(theta), cos(theta)]])
T_united[3:5, 3:5] = Matrix([[cos(theta), sin(theta)], 
                             [-sin(theta), cos(theta)]])
T_united[2, 2] = 1
T_united[5, 5] = 1

print('>>> Transition Stffness Matrix:')
eq_disp('T_{united}', T_united)
print('~~~~~~~~~~~')
eq_disp('U', 'TU^{(e)}')
eq_disp('F', 'TF^{(e)}')
eq_disp('KU', 'F')
eq_disp('K^{(e)}U^{(e)}', 'F^{(e)}')
eq_disp('K^{(e)}T^{-1}U', 'T^{-1}F')
eq_disp('TK^{(e)}T^{-1}U', 'F')
eq_disp('TK^{(e)}T^{-1}', 'K')
eq_disp('K^{(e)}', 'T^{-1}KT')
print('~~~~~~~~~~~')
eq_disp('T\'U', 'U^{(e)}')
eq_disp('T\'F', 'F^{(e)}')
eq_disp('T\'', 'T^{-1}')
eq_disp('T\'^{-1}K^{(e)}T\'', 'K')
eq_disp('K^{(e)}', 'T\'KT\'^{-1}')

