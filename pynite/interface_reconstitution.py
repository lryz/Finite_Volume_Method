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
    argsleft=eq[0]
    argsright=eq[1]
    coefleft=[]
    for i in range(0, len(argsleft)):
        coefleft.append([])
        for j in range(0,len(argsleft[i])):
            if type(argsleft[i][j])!=sp.core.function.Derivative :
                coefleft[i].append(argsleft[i][j])
    Flux1=[]
    Flux2=[]
    for i in range(0, len(argsright)):
        Flux1.append([])
        Flux2.append([])
        for j in range(0,len(argsright[i])):
            if type(argsright[i][j])==sp.core.function.Derivative :
                if argsright[i][j].args[1][0] == x :
                    Flux1[i].append((argsright[i][j].args[0].subs(
                        {xiplusundemi : xiplusun}) 
                        - argsright[i][j].args[0].subs(
                            {xiplusundemi : xi}))/(xiplusun-xi))
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
            Flux1[i][0]= Flux1[i][0]/2
            Flux2[i][0]= Flux2[i][0]/2
        Flux.append(Flux1[i])
        Flux.append(Flux2[i])
    out = recomposition([argsleft,Flux])
    return out