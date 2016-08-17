
# coding: utf-8

# In[1]:

from sympy import *
init_session(quiet=True, use_latex='mathjax', use_unicode=True)
from IPython.display import HTML, display

import os, sys
lib_path = os.path.abspath(os.path.join('..', 'src'))
sys.path.append(lib_path)

from EqDisp import eq_disp, eqs_disp


# <H4>The Wikipedia: <a href="https://en.wikipedia.org/wiki/Gaussian_quadrature" title="Gaussian Quadrature">click</a>

# In numerical analysis, a quadrature rule is an approximation of the definite integral of a function, usually stated as a weighted sum of function values at specified points within the domain of integration. (See numerical integration for more on quadrature rules.) An n-point <B>Gaussian quadrature rule</B>, named after Carl Friedrich Gauss, is a quadrature rule constructed to yield an exact result for polynomials of degree 2n − 1 or less by a suitable choice of the points xi and weights wi for i = 1, ..., n. The domain of integration for such a rule is conventionally taken as [−1, 1], so the rule is stated as
# 
# $${\int _{{-1}}^{1}f(x)\,dx=\sum _{{i=1}}^{n}w_{i}f(x_{i}).}$$
# 
# Gaussian quadrature as above will only produce accurate results if the function ${f(x)}$ is well approximated by a polynomial function within the range [−1, 1]. The method is not, for example, suitable for functions with singularities. However, if the integrated function can be written as ${f(x)=\omega (x)g(x)}$, where ${g(x)}$ is approximately polynomial and ${\omega (x)}$ is known, then alternative weights ${w_{i}'}$ and points ${x_{i}'}$ that depend on the weighting function ${\omega (x)}$ may give better results, where
# 
# $${\int _{{-1}}^{1}f(x)\,dx=\int _{{-1}}^{1}\omega (x)g(x)\,dx\approx \sum _{{i=1}}^{n}w_{i}'g(x_{i}').}$$
# 
# Common weighting functions include ${\omega (x)=1/{\sqrt {1-x^{2}}}}$ (Chebyshev–Gauss) and ${\omega (x)=e^{{-x^{2}}}}$ (Gauss–Hermite).
# 
# It can be shown (see Press, et al., or Stoer and Bulirsch) that the evaluation points xi are just the roots of a polynomial belonging to a class of orthogonal polynomials.

# <H4>The Book:

# $${I=\int _{{a}}^{b}f(x)\,dx=\sum _{{i=1}}^{n}w_{i}f(x_{i}).}$$
# 
# $${x=c_0+c_1\lambda}$$
# 
# $${a=c_0+c_1(-1)}$$
# $${a=c_0+c_1(1)}$$
# 
# $${c_0=\frac {b+a}{2}}$$
# $${c_1=\frac {b-a}{2}}$$
# 
# $${x=\frac {b+a}{2} + \frac {b-a}{2}\lambda}$$
# 
# $${dx=\frac {b-a}{2}d\lambda}$$
# 
# $${I=\int _{{-1}}^{1}f(\lambda)\,d\lambda=\sum _{{i=1}}^{n}w_{i}f(\lambda_{i}).}$$

# In[2]:

var('w_1 w_2 w_3 w_4 lambda_1 lambda_2 lambda_3 lambda_4')
_lambda_ = symbols('lambda');


# In[3]:

Eq_1_r = Integral(_lambda_**0, [_lambda_, -1, 1]); Eq_1_r_ = Eq_1_r.doit()
Eq_2_r = Integral(_lambda_**1, [_lambda_, -1, 1]); Eq_2_r_ = Eq_2_r.doit()
Eq_3_r = Integral(_lambda_**2, [_lambda_, -1, 1]); Eq_3_r_ = Eq_3_r.doit()
Eq_4_r = Integral(_lambda_**3, [_lambda_, -1, 1]); Eq_4_r_ = Eq_4_r.doit()
Eq_5_r = Integral(_lambda_**4, [_lambda_, -1, 1]); Eq_5_r_ = Eq_5_r.doit()
Eq_6_r = Integral(_lambda_**5, [_lambda_, -1, 1]); Eq_6_r_ = Eq_6_r.doit()
Eq_7_r = Integral(_lambda_**6, [_lambda_, -1, 1]); Eq_7_r_ = Eq_7_r.doit()
Eq_8_r = Integral(_lambda_**7, [_lambda_, -1, 1]); Eq_8_r_ = Eq_8_r.doit()


# <H5>1-Point Rule:

# In[7]:

Eq_1_l = w_1 * lambda_1**0
Eq_2_l = w_1 * lambda_1**1

eqs_disp([Eq_1_l, Eq_1_r, Eq_1_r_],
         [Eq_2_l, Eq_2_r, Eq_2_r_], inline=False)

sol = solve([Eq_1_l - Eq_1_r_, 
             Eq_2_l - Eq_2_r_])

w_1_val, lambda_1_val = sol[0][w_1], sol[0][lambda_1]

eqs_disp([w_1, w_1_val], 
         [lambda_1, lambda_1_val])


# <H5>2-Point Rule:

# In[4]:

Eq_1_l = w_1 * lambda_1**0 + w_2 * lambda_2**0
Eq_2_l = w_1 * lambda_1**1 + w_2 * lambda_2**1
Eq_3_l = w_1 * lambda_1**2 + w_2 * lambda_2**2
Eq_4_l = w_1 * lambda_1**3 + w_2 * lambda_2**3

eqs_disp([Eq_1_l, Eq_1_r, Eq_1_r_],
         [Eq_2_l, Eq_2_r, Eq_2_r_],
         [Eq_3_l, Eq_3_r, Eq_3_r_],
         [Eq_4_l, Eq_4_r, Eq_4_r_], inline=False)

sol = solve([Eq_1_l - Eq_1_r_, 
             Eq_2_l - Eq_2_r_, 
             Eq_3_l - Eq_3_r_, 
             Eq_4_l - Eq_4_r_])

w_1_val, w_2_val, lambda_1_val, lambda_2_val = sol[0][w_1], sol[0][w_2], sol[0][lambda_1], sol[0][lambda_2]

eqs_disp([w_1, w_1_val], 
         [w_2, w_2_val], 
         [lambda_1, lambda_1_val], 
         [lambda_2, lambda_2_val])


# <H5>3-Point Rule:

# In[5]:

Eq_1_l = w_1 * lambda_1**0 + w_2 * lambda_2**0 + w_3 * lambda_3**0
Eq_2_l = w_1 * lambda_1**1 + w_2 * lambda_2**1 + w_3 * lambda_3**1
Eq_3_l = w_1 * lambda_1**2 + w_2 * lambda_2**2 + w_3 * lambda_3**2
Eq_4_l = w_1 * lambda_1**3 + w_2 * lambda_2**3 + w_3 * lambda_3**3
Eq_5_l = w_1 * lambda_1**4 + w_2 * lambda_2**4 + w_3 * lambda_3**4
Eq_6_l = w_1 * lambda_1**5 + w_2 * lambda_2**5 + w_3 * lambda_3**5

eqs_disp([Eq_1_l, Eq_1_r, Eq_1_r_],
         [Eq_2_l, Eq_2_r, Eq_2_r_],
         [Eq_3_l, Eq_3_r, Eq_3_r_],
         [Eq_4_l, Eq_4_r, Eq_4_r_], 
         [Eq_5_l, Eq_5_r, Eq_5_r_], 
         [Eq_6_l, Eq_6_r, Eq_6_r_], inline=False)

sol = solve([Eq_1_l - Eq_1_r_, 
             Eq_2_l - Eq_2_r_, 
             Eq_3_l - Eq_3_r_, 
             Eq_4_l - Eq_4_r_,
             Eq_5_l - Eq_5_r_,
             Eq_6_l - Eq_6_r_])

w_1_val, w_2_val, w_3_val, lambda_1_val, lambda_2_val, lambda_3_val = sol[0][w_1], sol[0][w_2], sol[0][w_3], sol[0][lambda_1], sol[0][lambda_2], sol[0][lambda_3]

eqs_disp([w_1, w_1_val], 
         [w_2, w_2_val], 
         [w_3, w_3_val],
         [lambda_1, lambda_1_val], 
         [lambda_2, lambda_2_val],
         [lambda_3, lambda_3_val])


# <H5>4-Point Rule:

# In[6]:

Eq_1_l = w_1 * lambda_1**0 + w_2 * lambda_2**0 + w_3 * lambda_3**0 + w_4 * lambda_4**0
Eq_2_l = w_1 * lambda_1**1 + w_2 * lambda_2**1 + w_3 * lambda_3**1 + w_4 * lambda_4**1
Eq_3_l = w_1 * lambda_1**2 + w_2 * lambda_2**2 + w_3 * lambda_3**2 + w_4 * lambda_4**2
Eq_4_l = w_1 * lambda_1**3 + w_2 * lambda_2**3 + w_3 * lambda_3**3 + w_4 * lambda_4**3
Eq_5_l = w_1 * lambda_1**4 + w_2 * lambda_2**4 + w_3 * lambda_3**4 + w_4 * lambda_4**4
Eq_6_l = w_1 * lambda_1**5 + w_2 * lambda_2**5 + w_3 * lambda_3**5 + w_4 * lambda_4**5
Eq_7_l = w_1 * lambda_1**6 + w_2 * lambda_2**6 + w_3 * lambda_3**6 + w_4 * lambda_4**6
Eq_8_l = w_1 * lambda_1**7 + w_2 * lambda_2**7 + w_3 * lambda_3**7 + w_4 * lambda_4**7

eqs_disp([Eq_1_l, Eq_1_r, Eq_1_r_],
         [Eq_2_l, Eq_2_r, Eq_2_r_],
         [Eq_3_l, Eq_3_r, Eq_3_r_],
         [Eq_4_l, Eq_4_r, Eq_4_r_], 
         [Eq_5_l, Eq_5_r, Eq_5_r_], 
         [Eq_6_l, Eq_6_r, Eq_6_r_], 
         [Eq_7_l, Eq_7_r, Eq_7_r_], 
         [Eq_8_l, Eq_8_r, Eq_8_r_], inline=False)

sol = solve([Eq_1_l - Eq_1_r_, 
             Eq_2_l - Eq_2_r_, 
             Eq_3_l - Eq_3_r_, 
             Eq_4_l - Eq_4_r_,
             Eq_5_l - Eq_5_r_,
             Eq_6_l - Eq_6_r_,
             Eq_7_l - Eq_7_r_,
             Eq_8_l - Eq_8_r_])

w_1_val, w_2_val, w_3_val, w_4_val, lambda_1_val, lambda_2_val, lambda_3_val, lambda_4_val = sol[0][w_1], sol[0][w_2], sol[0][w_3], sol[0][w_4], sol[0][lambda_1], sol[0][lambda_2], sol[0][lambda_3], sol[0][lambda_4]

eqs_disp([w_1, w_1_val], 
         [w_2, w_2_val], 
         [w_3, w_3_val],
         [w_4, w_4_val],
         [lambda_1, lambda_1_val], 
         [lambda_2, lambda_2_val],
         [lambda_3, lambda_3_val],
         [lambda_4, lambda_4_val])


# In[ ]:



