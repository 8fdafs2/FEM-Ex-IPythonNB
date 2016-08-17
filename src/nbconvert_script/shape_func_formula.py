
# coding: utf-8

# In[1]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <H4>The Wikipedia: <a href="http://en.wikipedia.org/wiki/Lagrange_polynomial" title="Lagrange Polynomial">click</a>

# In numerical analysis, <B>Lagrange polynomials</B> are used for polynomial interpolation. For a given set of distinct points ${x_j}$ and numbers ${y_j}$, the <B>Lagrange polynomial</B> is the polynomial of the least degree that at each point ${x_j}$ assumes the corresponding value ${y_j}$ (i.e. the functions coincide at each point). The interpolating polynomial of the least degree is unique, however, and it is therefore more appropriate to speak of "the Lagrange form" of that unique polynomial rather than "the Lagrange interpolation polynomial", since the same polynomial can be arrived at through multiple methods. Although named after Joseph Louis Lagrange, who published it in 1795, it was first discovered in 1779 by Edward Waring and it is also an easy consequence of a formula published in 1783 by Leonhard Euler.
# 
# Lagrange interpolation is susceptible to Runge's phenomenon, and the fact that changing the interpolation points requires recalculating the entire interpolant can make Newton polynomials easier to use. <B>Lagrange polynomials</B> are used in the Newton–Cotes method of numerical integration and in Shamir's secret sharing scheme in cryptography.

# Given a set of k+1 data points
# 
# $${(x_{0},y_{0}),\ldots ,(x_{j},y_{j}),\ldots ,(x_{k},y_{k})}$$
# 
# where no two ${x_j}$ are the same, the interpolation polynomial in the Lagrange form is a linear combination
# 
# $${L(x):=\sum _{{j=0}}^{{k}}y_{j}\ell _{j}(x)}$$
# 
# of Lagrange basis polynomials
# 
# $${\ell _{j}(x):=\prod _{{{\begin{smallmatrix}0\leq m\leq k\\m\neq j\end{smallmatrix}}}}{\frac {x-x_{m}}{x_{j}-x_{m}}}={\frac {(x-x_{0})}{(x_{j}-x_{0})}}\cdots {\frac {(x-x_{{j-1}})}{(x_{j}-x_{{j-1}})}}{\frac {(x-x_{{j+1}})}{(x_{j}-x_{{j+1}})}}\cdots {\frac {(x-x_{k})}{(x_{j}-x_{k})}},}$$
# 
# where ${0\leq j\leq k}$. Note how, given the initial assumption that no two ${x_{i}}$ are the same, ${x_{j}-x_{m}\neq 0}$, so this expression is always well-defined. The reason pairs ${x_{i}=x_{j}}$ with ${y_{i}\neq y_{j}}$ are not allowed is that no interpolation function ${L}$ such that ${y_{i}=L(x_{i})}$ would exist; a function can only get one value for each argument ${x_i}$. On the other hand, if also ${y_{i}=y_{j}}$, then those two points would actually be one single point.
# 
# For all ${i\neq j}$, ${\ell _{j}(x)}$ includes the term ${(x-x_{i})}$ in the numerator, so the whole product will be zero at ${x=x_{i}}$:
# 
# $${\ell _{{j\neq i}}(x_{i})=\prod _{{m\neq j}}{\frac {x_{i}-x_{m}}{x_{j}-x_{m}}}={\frac {(x_{i}-x_{0})}{(x_{j}-x_{0})}}\cdots {\frac {(x_{i}-x_{i})}{(x_{j}-x_{i})}}\cdots {\frac {(x_{i}-x_{k})}{(x_{j}-x_{k})}}=0.}$$
# 
# On the other hand,
# 
# $${\ell _{i}(x_{i}):=\prod _{{m\neq i}}{\frac {x_{i}-x_{m}}{x_{i}-x_{m}}}=1}$$
# 
# In other words, all basis polynomials are zero at ${x=x_{i}}$, except ${\ell _{i}(x)}$, for which it holds that ${\ell _{i}(x_{i})=1}$, because it lacks the ${(x-x_{i})}$ term.
# 
# It follows that ${y_{i}\ell _{i}(x_{i})=y_{i}}$, so at each point ${x_{i}}$, ${L(x_{i})=y_{i}+0+0+\dots +0=y_{i}}$, showing that ${L}$ interpolates the function exactly.

# <H4>The Book:

# $${S_K = \prod _{M=1}^{N}\frac {X-X_M \,\,\, omit \,\,\, (X-X_K)}{X_K-X_M \,\,\, omit \,\,\, (X_K-X_K)}}$$

# In[2]:

var('X X_1 X_2 X_3 X_4 ell xi x');


# <H5>Linear Interpolation:

# In[3]:

sub_pairs_ = [(X_2 - X_1, ell)]
sub_pairs__ = [(X, x), (X_1, 0), (X_2, ell)]
sub_pairs___ = [(x, (xi + 1) * ell/2)]

S_1 = (X - X_2) / ((X_1 - X_2))
S_1_ = S_1.subs(sub_pairs_)
S_1__ = S_1_.subs(sub_pairs__).factor()
S_1___ = S_1__.subs(sub_pairs___).factor()

S_2 = (X - X_1) / ((X_2 - X_1))
S_2_ = S_2.subs(sub_pairs_)
S_2__ = S_2_.subs(sub_pairs__).factor()
S_2___ = S_2__.subs(sub_pairs___).factor()

eq_disp('S_1', S_1, S_1_, S_1__, S_1___)
eq_disp('S_2', S_2, S_2_, S_2__, S_2___)


# <H5>Quadratic Interpolation:

# In[4]:

sub_pairs_ = [(X_2 - X_1, ell/2), (X_3 - X_1, ell), (X_3 - X_2, ell/2)]
sub_pairs__ = [(X, x), (X_1, 0), (X_2, ell/2), (X_3, ell)]
sub_pairs___ = [(x, (xi + 1) * ell/2)]

S_1 = (X - X_2) * (X - X_3) / ((X_1 - X_2) * (X_1 - X_3))
S_1_ = S_1.subs(sub_pairs_)
S_1__ = S_1_.subs(sub_pairs__).factor()
S_1___ = S_1__.subs(sub_pairs___).factor()

S_2 = (X - X_1) * (X - X_3) / ((X_2 - X_1) * (X_2 - X_3))
S_2_ = S_2.subs(sub_pairs_)
S_2__ = S_2_.subs(sub_pairs__).factor()
S_2___ = S_2__.subs(sub_pairs___).factor()

S_3 = (X - X_1) * (X - X_2) / ((X_3 - X_1) * (X_3 - X_2))
S_3_ = S_3.subs(sub_pairs_)
S_3__ = S_3_.subs(sub_pairs__).factor()
S_3___ = S_3__.subs(sub_pairs___).factor()

eq_disp('S_1', S_1, S_1_, S_1__, S_1___)
eq_disp('S_2', S_2, S_2_, S_2__, S_2___)
eq_disp('S_3', S_3, S_3_, S_3__, S_3___)


# <H5>Cubic Interpolation:

# In[5]:

sub_pairs_ = [(X_2 - X_1, ell/3), (X_3 - X_1, ell/3*2), (X_3 - X_2, ell/3), (X_4 - X_1, ell), (X_4 - X_2, ell/3*2), (X_4 - X_3, ell/3)]
sub_pairs__ = [(X, x), (X_1, 0), (X_2, ell/3), (X_3, ell/3*2), (X_4, ell)]
sub_pairs___ = [(x, (xi + 1) * ell/2)]

S_1 = (X - X_2) * (X - X_3) * (X - X_4) / ((X_1 - X_2) * (X_1 - X_3) * (X_1 - X_4))
S_1_ = S_1.subs(sub_pairs_)
S_1__ = S_1_.subs(sub_pairs__).factor()
S_1___ = S_1__.subs(sub_pairs___).factor()

S_2 = (X - X_1) * (X - X_3) * (X - X_4) / ((X_2 - X_1) * (X_2 - X_3) * (X_2 - X_4))
S_2_ = S_2.subs(sub_pairs_)
S_2__ = S_2_.subs(sub_pairs__).factor()
S_2___ = S_2__.subs(sub_pairs___).factor()

S_3 = (X - X_1) * (X - X_2) * (X - X_4) / ((X_3 - X_1) * (X_3 - X_2) * (X_3 - X_4))
S_3_ = S_3.subs(sub_pairs_)
S_3__ = S_3_.subs(sub_pairs__).factor()
S_3___ = S_3__.subs(sub_pairs___).factor()

S_4 = (X - X_1) * (X - X_2) * (X - X_3) / ((X_4 - X_1) * (X_4 - X_2) * (X_4 - X_3))
S_4_ = S_4.subs(sub_pairs_)
S_4__ = S_4_.subs(sub_pairs__).factor()
S_4___ = S_4__.subs(sub_pairs___).factor()

eq_disp('S_1', S_1, S_1_, True, S_1__, S_1___)
eq_disp('S_2', S_2, S_2_, True, S_2__, S_2___)
eq_disp('S_3', S_3, S_3_, True, S_3__, S_3___)
eq_disp('S_4', S_4, S_4_, True, S_4__, S_4___)

