from IPython.display import HTML, display
from sympy import latex

#-----------------------------------------------------------------------------
# Basic Display Function
#-----------------------------------------------------------------------------


def htmlexpr_disp(htmlexpr):
    '''
    Display a HTML object
    '''
    display(htmlexpr)


def latexexpr_disp(latexexpr):
    '''
    Display a HTML object
    converted from a LaTeX equivalent presentation
    '''
    htmlexpr = HTML(latexexpr)
    display(htmlexpr)


def expr_disp(expr):
    '''
    Display a HTML object
    converted from a LaTeX equivalent presentation of
    a plain string
    '''
    latexexpr = latex(expr)
    htmlexpr = HTML(latexexpr)
    display(htmlexpr)

#-----------------------------------------------------------------------------
# Advanced Display Function
#-----------------------------------------------------------------------------


def eq_disp(*exprs, inline=False):
    '''
    Display one equation in inline/separate mode
    For each expression in expressions (exprs),
    obtain its LaTeX equivalent presentation first
    based on the inline/separate mode,
    then connect them via '=' for the display
    '''
    latexexpr = latexeq_get(*exprs, inline=inline)

    latexexpr_disp(latexexpr)


def eqs_disp(*exprspairs, inline=True, spacing=5):
    '''
    Display multiple equations in inline/separate mode
    For each expressions (exprs) in expressions-pairs (exprspairs),
        For each expression in expressions (exprs),
        obtain its LaTeX equivalent presentation first
        based on the inline/separate mode,
    connect them via inline mark or empty string for the display
    '''
    s_list = []

    for exprs in exprspairs:
        s_list.append(latexeq_get(*exprs, inline=inline))

    if inline:
        latexexpr = ('$' + r'\,' * spacing + '$').join(s_list)
    else:
        latexexpr = ''.join(s_list)

    latexexpr_disp(latexexpr)

#-----------------------------------------------------------------------------
# Conversion Function
#-----------------------------------------------------------------------------


def htmleq_get(*exprs, inline=False):
    '''
    Get HTML object of one equation in inline/separate mode
    For each expression in expressions (exprs),
    obtain its LaTeX equivalent presentation first
    based on the inline/separate mode,
    then connect them via '=',
    finally convert the result into a HTML object to return
    '''
    latexexpr = latexeq_get(*exprs, inline=inline)

    htmlexpr = HTML(latexexpr)

    return htmlexpr


def latexeq_get(*exprs, inline=False):
    '''
    Get LaTeX equivalent presentation of one equation in inline/separate mode
    For each expression in expressions (exprs),
    obtain its LaTeX equivalent presentation first
    based on the inline/separate mode,
    then connect them via '=' for the final presentation to return
    '''

    if exprs:
        exprs = list(exprs)

        while exprs:
            expr = exprs[0]
            if isinstance(expr, bool):
                exprs.pop(0)
            else:
                break

        if exprs:
            i = 1
            connectors = []
            while i <= len(exprs) - 1:
                expr = exprs[i]
                if isinstance(expr, bool) and expr is True:
                    connectors.append(r'\\&=')
                    exprs.pop(i)
                else:
                    connectors.append('=')
                i += 1

            i = 1
            while i <= len(connectors) - 1:
                if len(connectors[i]) == 4:
                    if len(connectors[i - 1]) == 1:
                        connectors[i - 1] = '&='
                i += 1

#             if len(connectors) == 1:
#                 connectors[0] = '='
#             else:
#                 i = 1
#                 while i <= len(connectors) - 1:
#                     if len(connectors[i]) == 4:
#                         if len(connectors[i - 1]) == 1:
#                             connectors[i - 1] = '&='
#                     i += 1

            s = '{' + latex(exprs[0]) + '}'
            for pair in zip(connectors, exprs[1:]):
                s += pair[0] + '{' + latex(pair[1]) + '}'

            if inline:
                s = r'\(' + s + r'\)'
            else:
                s = r'\begin{equation}\begin{aligned}' + \
                    s + '\end{aligned}\end{equation}'
#             print(connectors)
#             print(s)
            return s

    return
