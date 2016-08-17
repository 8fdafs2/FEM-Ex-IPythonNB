
# coding: utf-8

# In[65]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <img src="image/img_01.png" alt="img_01" align=left>

# For any element:
# 
# $${S_1 = -\frac {1}{\ell}(-\ell+x)}$$
# $${S_2 = \frac {x}{\ell}}$$

# In[66]:

var('ell x xi X X_1 X_2 S_1 S_2 T_1 T_2');


# In[67]:

T = S_1 * T_1 + S_2 * T_2

S_1_ = (X - X_2) / (X_1 - X_2)
S_2_ = (X - X_1) / (X_2 - X_1)

T_ = S_1_ * T_1 + S_2_ * T_2

eq_disp('T', T, T_)


# In[68]:

S_1__ = 1/ell*(ell-x)
S_2__ = x/ell

T__ = S_1__ * T_1 + S_2__ * T_2 

eq_disp('T', T__)


# ${X = 8 \,\,\, >>> \,\,\, \ell = 5 \,\,\, x = 3 \,\,\, T_1 = 34 \,\,\, T_2 = 20}$

# In[69]:

ell_val, x_val, T_1_val, T_2_val = 5, 3, 34, 20


# In[70]:

T_val = T__.subs([(ell, ell_val), (x, x_val), (T_1, T_1_val), (T_2, T_2_val)])

eq_disp('T', T_val.evalf())


# In[71]:

T___ = T__.subs([(x, (xi + 1) * ell/2)])

eq_disp('T', T___.simplify().factor(T_1, T_2))


# ${X = 7.5 \,\,\, >>> \,\,\, \xi = 2*x/\ell-1 = 2*2.5/5-1 = 0 \,\,\, T_1 = 34 \,\,\, T_2 = 20}$

# In[72]:

xi_val, T_1_val, T_2_val = 0, 34, 20


# In[73]:

T_val = T___.subs([(xi, xi_val), (T_1, T_1_val), (T_2, T_2_val)])

eq_disp('T', T_val.evalf())

