
# coding: utf-8

# In[99]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <H4>Integral Calculation:
# 
# $${I = \int_{2}^{6} (x^2+5x+3)dx}$$

# In[100]:

var('x f xi ell')
_lambda_ = symbols('lambda');


# Analytical Solution:

# In[101]:

f = x**2+5*x+3
I = Integral(f, [x, 2, 6])
eq_disp('I', I, I.doit(), I.doit().evalf())


# $${x=\frac{6+2}{2}+\frac{6-2}{2}\lambda=4+2\lambda}$$
# $${dx = 2d\lambda}$$

# In[102]:

f_ = f.subs(x, 4+2*_lambda_) * 2
I_ = Integral(f_, [_lambda_, -1, 1])
eq_disp('I', I_, I_.doit(), I_.doit().evalf())


# In[103]:

I__ = 2 * f_.subs(_lambda_, 0)
I___ = 1 * f_.subs(_lambda_, -sqrt(3)/3) + 1 * f_.subs(_lambda_, sqrt(3)/3)
I____ = 8/9 * f_.subs(_lambda_, 0) + 5/9 * f_.subs(_lambda_, -sqrt(15)/5) + 5/9 * f_.subs(_lambda_, sqrt(15)/5)

print('1-Point:')
eq_disp('I', I__.evalf())
print('2-Point:')
eq_disp('I', I___.evalf())
print('3-Point:')
eq_disp('I', I____.evalf())


# <H4>Integral Calculation:
# 
# $${
# I = 
# \int_{X_i}^{X_j} S_j^2 dX = 
# \int_{X_i}^{X_j} {\left(\frac {X-X_i}{\ell}\right)}^2 dX = 
# \int_{0}^{\ell} {\left(\frac {x}{\ell}\right)}^2 dx = 
# \frac {\ell}{2} \int_{-1}^{1} {\left[\frac {1}{2}(1+\xi)\right]}^2 d\xi
# }$$

# Analytical Solution:

# In[104]:

f = ell / 2 * ((1 + xi) / 2)**2
I = Integral(f, [xi, -1, 1])
eq_disp('I', I, I.doit(), I.doit().evalf())


# In[105]:

f_ = f.subs(xi, _lambda_)
I_ = Integral(f_, [_lambda_, -1, 1])
eq_disp('I', I_, I_.doit(), I_.doit().evalf())


# In[106]:

I__ = 2 * f_.subs(_lambda_, 0)
I___ = 1 * f_.subs(_lambda_, -sqrt(3)/3) + 1 * f_.subs(_lambda_, sqrt(3)/3)
I____ = 8/9 * f_.subs(_lambda_, 0) + 5/9 * f_.subs(_lambda_, -sqrt(15)/5) + 5/9 * f_.subs(_lambda_, sqrt(15)/5)

print('1-Point:')
eq_disp('I', I__.evalf())
print('2-Point:')
eq_disp('I', I___.evalf())
print('3-Point:')
eq_disp('I', I____.evalf())


# In[ ]:



