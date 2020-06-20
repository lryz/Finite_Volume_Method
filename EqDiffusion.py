import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from random import *
import numpy as np
from sympy.interactive import init_printing
from sympy.interactive import printing
x = sp.Symbol('x')
t = sp.Symbol('t')
u = sp.Function('u')
a = sp.Symbol('x_i-1/2', constant = True)
b = sp.Symbol('x_i+1/2')
c = sp.Symbol('c')

def equation_diffusion(D):
    """
    La fonction equation_diffusion permet de générer un objet 
    de classe <class 'sympy.core.relational.Equality'>
    qui représente l'équation de diffusion.

    Elle prend en entrée le coefficient de diffusion D.
    """
    u = sp.Function('u')
    x = sp.Symbol('x')
    t = sp.Symbol('t')
    eq1 = sp.Eq(u(x,t).diff(t), D*u(x,t).diff(x).diff(x))
    return eq1


def maillage1D_regulier(debut,fin,pas):
    N=int((fin-debut)/pas+1)
    milieu=np.linspace(debut,fin,N)
    return milieu

def initialisation_aleatoire(nbr_pts):
    u=np.ones(nbr_pts)
    for i in range(0,nbr_pts):
        u[i]=u[i]+random()
    return u 

def pre(expr,ls):
    """
    Cette fonction nous permet de ressortir la liste des arguments 'importants' de notre equation.
    On va avoir notamment une liste contenant la valeur des différents coefficients réels
    ainsi que les différentes dérivées partielles.
    Cette fonction est récursive puisque la représentation des equations de sympy se fait selon un arbre binaire 
    (voir https://docs.sympy.org/latest/tutorial/manipulation.html).
    Les conditions nous permettent de nous arreter aux valeurs importantes que l'on veut garder et ainsi éviter de les redécomposer.
    Notamment dans le cas des dérivées qui contiennent différents arguments.
    #La condition est probablement inutile pour le float
    Ajustement possibles en rajoutant d'autre types d'arguments comme des fonctions. (notamment le u(1-u) utilisé dans Keller-Segel)
    """
    if (type(expr)==type(sp.Float(0.0))) :
        ls.append(expr)
    elif (type(expr)==type(u(x,t).diff(x))):
        ls.append(expr)
    else :
        for arg in expr.args:
            pre(arg,ls)
            #On rappelle la fonction pour décomposer encore plus nos arguments.

def flux_differences_finies(equation):
    """
    La fonction permet d'intégrer sur x l'équation de diffusion de manière "Formelle".
    On va séparer les arguments de chaque côté de l'équation en prenant les listes des arguments.
    On est obligé de créer des listes bis car on ne peut pas assigner des valeurs à des variables de type "Derivative".

    #WIP#
    Pour l'instant, on ne peut pas rajouter de terme d'addition dans l'équation. On est obligé d'avoir des produits de chaque cotés.
    """
    eqleft = equation.lhs
    eqright = equation.rhs
    listleft=[]
    listright=[]
    pre(eqleft,listleft)
    pre(eqright, listright)
    listleftbis=[]
    for i in range(0,len(listleft)) :
        if (type(listleft[i])==type(u(x,t).diff(x))):
            if listleft[i].args[1][0] == t :
                listleftbis.append((b-a)*listleft[i].args[0].diff(listleft[i].args[1]))
            if listleft[i].args[1][0] == x :
                listleftbis.append(sp.Derivative(listleft[i].args[0].subs({x : b})
                    ,(listleft[i].args[1])[0],(listleft[i].args[1])[1]-1)
                    - sp.Derivative(listleft[i].args[0].subs({x : a}),
                        (listleft[i].args[1])[0],(listleft[i].args[1])[1]-1))
        elif (type(listleft[i])==type(sp.Float(0.0))) :
            listleftbis.append(listleft[i])
        else :
            listleftbis.append(listleft[i])

    listrightbis=[]

    for i in range(0,len(listright)) :
        if (type(listright[i])==type(u(x,t).diff(x))):
            if listright[i].args[1][0] == t :
                listrightbis.append((b-a)*listright[i].args[0].diff(listright[i].args[1]))
            if listright[i].args[1][0] == x :
                listrightbis.append(sp.Derivative(listright[i].args[0].subs({x : b})
                    ,(listright[i].args[1])[0],(listright[i].args[1])[1]-1)
                    - sp.Derivative(listright[i].args[0].subs({x : a}),
                        (listright[i].args[1])[0],(listright[i].args[1])[1]-1))
        elif (type(listright[i])==type(sp.Float(0.0))) :
            listrightbis.append(listright[i])
        else :
            listrightbis.append(listright[i])

    return sp.Eq(np.prod(listleftbis), np.prod(listrightbis))





eq = equation_diffusion(0.5)
sp.pprint(eq)
maillage = maillage1D_regulier(0,1,0.02)
uinit=initialisation_aleatoire(len(maillage))
test = sp.Eq(x*t+a, b + u(x,t).diff(x))
expr = flux_differences_finies(eq)
print(sp.printing.latex(expr))
#for arg in sp.preorder_traversal(test):
#    print(arg)
sp.preview(expr, viewer='file', filename='output.png')