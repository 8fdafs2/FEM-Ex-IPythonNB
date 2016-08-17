
# coding: utf-8

# In[12]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/4.7.png" alt="Figure 4.7" align=left>

# In[13]:

# Global
var('c_1 c_2 Y l y')
var('u_i u_j Y_i Y_j F_i F_j')
var('sigma epsilon V E A k');


# In[14]:

# Displacement Represented by Shape Function
print('>>> Displacement Represented by Shape Function:')
eq_disp('u_e', c_1 + c_2 * Y)
print('~~~~~~~~~~~')
eqs_disp(['u_i', c_1 + c_2 * Y_i], ['u_j', c_1 + c_2 * Y_j], inline=False)

Eq_1 = c_1 + c_2 * Y_i - u_i
Eq_2 = c_1 + c_2 * Y_j - u_j

sol = solve([Eq_1, Eq_2], [c_1, c_2])

c_1, c_2 = sol[c_1], sol[c_2]

print('~~~~~~~~~~~')
eqs_disp(['c_1', c_1], 
         ['c_2', c_2], 
         inline=False)

u_e = c_1 + c_2 * Y
eq_disp('u_e', u_e)

u_e = expand(u_e)

S_i = simplify(u_e.coeff(u_i))
S_j = simplify(u_e.coeff(u_j))

S_i_ = S_i.subs(Y_j - Y_i, l)
S_j_ = S_j.subs(Y_j - Y_i, l)

S_i__ = S_i_.subs(Y_j - Y, l - y)
S_j__ = S_j_.subs(Y - Y_i, y)

u_e_vec = Matrix([u_i, u_j])

print('~~~~~~~~~~~')
eq_disp('u_e',
        'S_iu_i + S_ju_j',
        MatMul(Matrix(1, 2, ['S_i', 'S_j']), u_e_vec),
        MatMul(Matrix(1, 2, [S_i, S_j]), u_e_vec), True,
        MatMul(Matrix(1, 2, [S_i_, S_j_]), u_e_vec), True,
        MatMul(Matrix(1, 2, [S_i__, S_j__]), u_e_vec))

S_i = S_i__
S_j = S_j__

u_e = simplify(S_i * u_i + S_j * u_j)

print('~~~~~~~~~~~')

eq_disp('u_e',
        'S_iu_i + S_ju_j',
        u_e)


# In[15]:

# Strain Energy
epsilon_ = Derivative('u_e', y)
epsilon__ = Derivative(u_e, y)
epsilon___ = epsilon__.doit()

Lambda = Integral(sigma * epsilon / 2, [V,])
Lambda_ = Lambda.subs(sigma, E * epsilon)
Lambda__ = Lambda_.subs(epsilon, epsilon___)
Lambda___ = E * epsilon___ ** 2 / 2 * A * l

print('>>> Strain Energy:')
eq_disp(r'\epsilon',
        epsilon_,
        epsilon__,
        epsilon___)
print('~~~~~~~~~~~')
eq_disp(r'\Lambda^{(e)}',
        Lambda,
        Lambda_,
        Lambda__,
        Lambda___)


# In[16]:

# Stiffness Matrix
print('>>> Stffness Matrix:')
eq_disp(r'\Pi',
        r'\sum_{e=1}^{n} \Lambda^{(e)} - \sum_{k=1}^{m} F_ku_k')

eq_disp(r'\frac{\partial \Pi}{\partial u_k}',
        r'\frac{\partial}{\partial u_k} \sum_{e=1}^{n} \Lambda^{(e)} - \frac{\partial}{\partial u_k} \sum_{i=k}^{m} F_ku_k',
        r'0, \,\,\,\,\,\, k = 1,2,3,...,m')

eq_disp(Matrix([Derivative(symbols(r'\Lambda'), u_i), Derivative(symbols(r'\Lambda'), u_j)]),
        Matrix([simplify(diff(Lambda___, u_i).doit()), simplify(diff(Lambda___, u_j).doit())]),
        MatMul(MatMul(A * E / l, Matrix([[1, -1], [-1, 1]])), Matrix([u_i, u_j])), True,
        MatMul(Matrix([[k, -k], [-k, k]]), Matrix([u_i, u_j])), True,
        'K^{(e)}U^{(e)}', True,
        Matrix([F_i, F_j]),
        'F^{(e)}')


# In[ ]:



