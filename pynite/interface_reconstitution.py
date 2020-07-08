import sympy as sp
from pynite.tools import *
x = sp.Symbol('x')
t = sp.Symbol('t')
u = sp.Function('u')
c = sp.Function('c')
phi = sp.Function('\phi')
ximoinsundemi = sp.Symbol('x_i-1/2')
xiplusundemi = sp.Symbol('x_i+1/2')
xi = sp.Symbol('x_i')
xiplusun = sp.Symbol('x_i+1')
ximoinsun = sp.Symbol('x_i-1')

def finite_difference(eq):
    """
    Interface reconstituion using finite difference method (first order)

    Parameters
    ----------
    eq: list 
          decomposed equation
        
    Returns
    -------
        decomposed equation with terms like d(u(xi+1/2))/dx replaced by (u(x+1)-u(xi))/(xi+1-xi)
    """
    argsleft=eq[0]
    argsright=eq[1]
    coefleft=[]
    for i in range(0, len(argsleft)):
        coefleft.append([])
        for j in range(0,len(argsleft[i])):
            if type(argsleft[i][j])!=sp.core.function.Derivative :
                #if our term isn't a derivative we don't modify it
                coefleft[i].append(argsleft[i][j])

    Flux1=[]    #Output flux (we consider xi+1/2)
    Flux2=[]    #Input flux (we condider xi-1/2)
    for i in range(0, len(argsright)):
        Flux1.append([])
        Flux2.append([])
        for j in range(0,len(argsright[i])):
            if type(argsright[i][j])==sp.core.function.Derivative :
                if argsright[i][j].args[1][0] == x :
                    #x derivative 
                    #d(u(xi+1/2))/dx replaced by (u(xi+1)-u(xi))/(xi+1-xi)
                    Flux1[i].append((argsright[i][j].args[0].subs(
                        {xiplusundemi : xiplusun}) 
                        - argsright[i][j].args[0].subs(
                            {xiplusundemi : xi}))/(xiplusun-xi))
                    #d(u(xi-1/2))/dx replaced by (u(xi)-u(xi-1))/(xi-xi-1)
                    Flux2[i].append(
                        (argsright[i][j].args[0].subs({ximoinsundemi : xi}) 
                        - argsright[i][j].args[0].subs(
                            {ximoinsundemi : ximoinsun}))/(xi-ximoinsun))
            else :
                Flux1[i].append(argsright[i][j])
                Flux2[i].append(argsright[i][j])
    Flux = []
    for i in range(0,len(Flux1)):
        if Flux1[i]==Flux2[i]:
            #in case we have constant terms in our addition we divide by 2 on each side to keep, the same coefficient after addition of each flux
            Flux1[i][0]= Flux1[i][0]/2
            Flux2[i][0]= Flux2[i][0]/2
        Flux.append(Flux1[i])
        Flux.append(Flux2[i])
    out = recomposition([argsleft,Flux])
    return out