
# coding: utf-8

# In[15]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/4.8.png" alt="Figure 4.8" align=left>
# <img src="image/4.9.png" alt="Figure 4.9" align=left>

# In[16]:

# Global
var('c_1 c_2 c_3 c_4 l x')
var('U_i1 U_i2 U_j1 U_j2')
var('w R_1 M_1 R_2 M_2')
var('sigma epsilon V E A I D U');


# In[17]:

# Displacement Represented by Shape Function
print('>>> Displacement Represented by Shape Function:')
eq_disp('v_e', c_1 + c_2 * x + c_3 * x ** 2 + c_4 * x ** 3)
print('~~~~~~~~~~~')
eqs_disp(['U_{i1}', c_1], 
         ['U_{i2}', c_2],
         ['U_{j1}', c_1 + c_2 * l + c_3 * l ** 2 + c_4 * l ** 3],
         ['U_{j2}', c_2 + 2 * c_3 * l + 3 * c_4 * l ** 2],
         inline=False)

Eq_1 = c_1 - U_i1
Eq_2 = c_2 - U_i2
Eq_3 = c_1 + c_2 * l + c_3 * l ** 2 + c_4 * l ** 3 - U_j1
Eq_4 = c_2 + 2 * c_3 * l + 3 * c_4 * l ** 2 - U_j2

sol = solve([Eq_1, Eq_2, Eq_3, Eq_4], [c_1, c_2, c_3, c_4])

c_1, c_2, c_3, c_4 = sol[c_1], sol[c_2], sol[c_3], sol[c_4]

print('~~~~~~~~~~~')
eqs_disp(['c_1', c_1], 
         ['c_2', c_2], 
         ['c_3', c_3],
         ['c_4', c_4],
         inline=False)

v_e = c_1 + c_2 * x + c_3 * x ** 2 + c_4 * x ** 3
eq_disp('v_e', v_e)

v_e = expand(v_e)

S_i1 = simplify(v_e.coeff(U_i1))
S_i2 = simplify(v_e.coeff(U_i2))
S_j1 = simplify(v_e.coeff(U_j1))
S_j2 = simplify(v_e.coeff(U_j2))

S_e_vec = Matrix(1, 4, [S_i1, S_i2, S_j1, S_j2])
v_e_vec = Matrix([U_i1, U_i2, U_j1, U_j2])

print('~~~~~~~~~~~')
eq_disp('v_e',
        'S_{i1}U_{i1} + S_{i2}U_{i2} + S_{j1}U_{j1} + S_{j2}U_{j2}', True,
        MatMul(Matrix(1, 4, [symbols('S_{i1}'), symbols('S_{i2}'), symbols('S_{j1}'), symbols('S_{j2}')]), v_e_vec), True,
        MatMul(S_e_vec, v_e_vec))


# In[18]:

# Strain Energy
epsilon_ = -y * Derivative('v_e', x, x)
epsilon__ = -y * Derivative(v_e, x, x)
epsilon___ = epsilon__.doit()
epsilon____ = expand(epsilon___)

D_i1 = simplify(epsilon____.coeff(U_i1)) / -y
D_i2 = simplify(epsilon____.coeff(U_i2)) / -y
D_j1 = simplify(epsilon____.coeff(U_j1)) / -y
D_j2 = simplify(epsilon____.coeff(U_j2)) / -y

D_e_vec = Matrix(1, 4, [D_i1, D_i2, D_j1, D_j2])

Lambda = Integral(sigma * epsilon / 2, [V,])
Lambda_ = Lambda.subs(sigma, E * epsilon)
Lambda__ = Lambda_.subs(epsilon, epsilon_)
Lambda___ = E / 2 * Integral((epsilon_ / y) ** 2, [x, 0, l]) * Integral(y ** 2, [A,])
Lambda____ = Lambda___.subs(Integral(y ** 2, [A,]), I)

print('>>> Strain Energy:')
eq_disp(r'\epsilon',
        epsilon_,
        epsilon__, True,
        epsilon___, True,
        MatMul(-y, D_e_vec, v_e_vec), True,
        MatMul(-y, Matrix(1, 4, [symbols('D_{i1}'), symbols('D_{i2}'), symbols('D_{j1}'), symbols('D_{j2}')]), v_e_vec), True,
        '-y\{D\}\{U\}'
        )
print('~~~~~~~~~~~')
eq_disp(r'\Lambda^{(e)}',
        Lambda,
        Lambda_, 
        Lambda__, True,
        Lambda___, True,
        Lambda____, True,
        r'\frac{EI}{2}\int_0^{l}(\{D\}\{U\})(\{D\}\{U\})dx', True,
        r'\frac{EI}{2}\int_0^{l}\{U\}^T\{D\}^T\{D\}\{U\}dx')


# In[19]:

# Stiffness Matrix
d_Lambda_d_U = MatMul(E * I, Integral(D_e_vec.transpose() * D_e_vec, [x, 0, l]), v_e_vec)
d_Lambda_d_U_ = MatMul(E * I/(l**3), (l**3) * integrate(D_e_vec.transpose() * D_e_vec, [x, 0, l]), v_e_vec)

print('Stffness Matrix:')
eq_disp(r'\Pi',
        r'\sum_{e=1}^{n} \Lambda^{(e)} - \sum_{k=1}^{m} F_kU_k')

eq_disp(r'\frac{\partial \Pi}{\partial U_k}',
        r'\frac{\partial}{\partial U_k} \sum_{e=1}^{n} \Lambda^{(e)} - \frac{\partial}{\partial U_k} \sum_{i=k}^{m} F_kU_k',
        r'0, \,\,\,\,\,\, k = 1,2,3,...,m')

eq_disp(r'\frac{\partial\Lambda^{(e)}}{\partial U^{(e)}}',
        r'EI\int_0^{l}\{D\}^T\{D\}dxU^{(e)}', True,
        d_Lambda_d_U, True,
        d_Lambda_d_U_,
        r'K^{(e)}U^{(e)}')


# In[22]:

# External Force Vector
Eq_1 = w * l**3 / 6 + R_1 * l**2 / 2 - M_1 * l
Eq_2 = w * l**4 / 24 + R_1 * l**3 / 6 - M_1 * l**2 / 2
sol = solve([Eq_1, Eq_2], [R_1, M_1])
R_1_val, M_1_val = sol[R_1], sol[M_1]
R_2_val, M_2_val = R_1_val, -M_1_val
F_e_vec = Matrix(4, 1, [-R_1_val, -M_1_val, -R_2_val, -M_2_val])

print('>>> External Force Vector:')
eq_disp(r'EI\frac{d^4 v}{dx^4}',
        r'\frac{dV(x)}{dx}',
        r'-w')
print('~~~~~~~~~~~')
eq_disp('V(0)', 'c_1', 'R_1')
eq_disp(r'EI\frac{d^3 v}{dx^3}',
        r'V(x)',
        r'wx+c_1',
        r'wx+R_1')
print('~~~~~~~~~~~')
eq_disp('M(0)', 'c_2', 'M_1')
eq_disp(r'EI\frac{d^2 v}{dx^2}',
        r'M(x)',
        r'wx^2/2+c_1x+c_2',
        r'wx^2/2+R_1x+M_1')
print('~~~~~~~~~~~')
eq_disp(r'\theta(0)', 'c_3', '0')
eq_disp(r'EI\frac{d v}{dx}',
        r'EI\theta(x)',
        r'wx^3/6+c_1x^2/2+c_2x+c_3',
        r'wx^3/6+R_1x^2/2+M_1x')
print('~~~~~~~~~~~')
eq_disp('v(0)', 'c_4', '0')
eq_disp(r'EIv',
        r'EIv(x)',
        r'wx^4/24+c_1x^3/6+c_2x^2/2+c_3x+c_4',
        r'wx^4/24+R_1x^3/6+M_1x^2/2')
print('~~~~~~~~~~~')
eq_disp(r'\theta(l)', '0')
eq_disp(r'EI\theta(l)',
        r'wl^3/6+R_1l^2/2+M_1l',
        0)
eq_disp('v(l)', '0')
eq_disp(r'EIv(l)',
        r'wl^4/24+R_1l^3/6+M_1l^2/2',
        0)
print('~~~~~~~~~~~')
eq_disp(Matrix(4, 1, ['-R_1', '-M_1', '-R_2', '-M_2']),
        F_e_vec,
        Matrix(4, 1, ['F_i1', 'F_i2', 'F_j1', 'F_j2']),
        '\{F\}^{(e)}')


# In[21]:

# Equilibrium Equation
print('>>> Equilibrium Equation:')
eq_disp(r'K^{(e)}U^{(e)}', 'F^{(e)}')
eq_disp(d_Lambda_d_U_, F_e_vec)

