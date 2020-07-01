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
a = sp.Symbol('x_i-1/2')
b = sp.Symbol('x_i+1/2')
c = sp.Symbol('x_i')
d = sp.Symbol('x_i+1')
e = sp.Symbol('x_i-1')
tn = sp.Symbol('t_n')
Deltat = sp.Symbol('\Delta t')
uni = sp.Symbol('u__n_i')
unplusuni = sp.Symbol('u_i__n+1')
uniplusun = sp.Symbol('u_i+1__n')
unplusuniplusun = sp.Symbol('u_i+1__n')
unimoinsun = sp.Symbol('u_i-1__n')
U = sp.Symbol('U')
Un = sp.Symbol('U__n')
Unplusun = sp.Symbol('U__n+1')

"""
Pour généraliser les fonctions on doit pouvoir rajouter une condition vérifiant si 
le type est une somme ou une multiplication et utiliser une décomposition puis reconstruction en conséquence.
Pour l'equation de diffusion comme on a que des produits (sans la source pour l'instant) on n'en a pas encore besoin 
"""

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
    maillage=np.linspace(debut,fin,N)
    return maillage

def initialisation_aleatoire(nbr_pts):
    u=np.ones(nbr_pts)
    for i in range(0,nbr_pts):
        u[i]=u[i]+random()
    return u eq1eq1

def calcul_milieu(debut,fin):
    return (debut+fin)/2

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
    elif (type(expr)==type(x)):
        ls.append(expr)
    elif (type(expr)==sp.core.numbers.NegativeOne):
        ls.append(expr)
        #Le problème ici c'est juste qu'on rajoute un -1 dans la liste. 
        #On a toujours le problème de savoir si c'est une addition ou une multiplication
        #Une fois que la liste est créeer on peut appeler la fonction simpl_moinsun(ls)
    else :
        for arg in expr.args:
            pre(arg,ls)
            #On rappelle la fonction pour décomposer encore plus nos arguments.


def simpl_moinsun(ls):
    """
    Cette fonction sert à reconstituer les valeurs négatives à partir des arguments
    """
    for i in range(0,len(ls)-1):
        if ls[i]==-1 :
            ls[i]=ls[i]*ls[i+1]
            del ls[i+1]
            

def flux_dimension1(equation):
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


def reconstitution_differences_finies(equation):
    """
    Calcul du flux à l'interface.
    """
    tmp = equation.rhs.args
    F = tmp[0]
    argsF = []
    pre(F,argsF)
    argsFluxnumerique1 = []
    argsFluxnumerique2 = []
    for i in range(0,len(argsF)):
        if type(argsF[i]) == type(u(x,t).diff(x)) :
            if argsF[i].args[1][0] == x :
                argsFluxnumerique1.append((argsF[i].args[0].subs({b : d}) - argsF[i].args[0].subs({b : c}))/(d-c))
                argsFluxnumerique2.append((argsF[i].args[0].subs({b : c}) - argsF[i].args[0].subs({b : e}))/(c-e))
        else :
            argsFluxnumerique1.append(argsF[i])
            argsFluxnumerique2.append(argsF[i])

    Flux1 = np.prod(argsFluxnumerique1)
    Flux2 = np.prod(argsFluxnumerique2)

    tmp = equation.lhs.args
    coef=[]
    left=[]
    for i in tmp :
        if type(i) != type(u(x,t).diff(x)) :
            coef.append(1/i)
        else :
            left.append(i)
    return sp.Eq(np.prod(left), np.prod(coef)*(Flux1-Flux2))

def schema_explicite(equation):
    """
    On suppose que l'équation est "triée" c'est à dire que la partie gauche contient les dérivées par rapport à t
    Cette fonction va renvoyer le schéma explicite sous forme "factorisée" 
    (c'est-à-dire qu'on aura : u^{n+1}_i = (facteur1 * u^n_{i+1}+ facteur2 * u^n_i + facteur3 * u^n_{i-1})/dividende)
    On va donc ensuite pouvoir la passer sous forme matricielle
    """
    left = 0
    if type(equation.lhs) == type(u(x,t).diff(x)) :
        if equation.lhs.args[1][0] == t :
            left = (unplusuni-uni)/Deltat
    rhs = equation.rhs.subs({t : tn})
    rhs = rhs.subs({u(c,tn) : uni})
    rhs = rhs.subs({u(d,tn) : uniplusun})
    rhs = rhs.subs({u(e,tn) : unimoinsun})
    argsleft = []
    for i in left.args :
        if type(i) == type(x**2):
            rhs = i**(-1)*rhs
        else :
            tmp = []
            pre(i,tmp)
            simpl_moinsun(tmp)
            for j in tmp:
                if j == unplusuni:
                    argsleft.append(j)
                else :
                    rhs = rhs - j
    return sp.Eq(np.prod(argsleft),rhs.factor(uniplusun,unimoinsun,uni))


def calcul_coeffs(equation):
    """
    Cette fonction prend en entrée notre équation factorisée après un schéma (uniquement Euler explicite pour l'instant) 
    et retourne une liste contenant les coefficients (dans l'ordre) des u^n_i-1, u^n_i, u^n_i+1
    """
    coefficients=[]
    for i in equation.rhs.args :
        if type(i)!= sp.core.add.Add :
            coefficients.append(i)
        elif type(i)==sp.core.add.Add :
            for tmp in i.args :
                if (type(tmp)==sp.core.mul.Mul):
                    if tmp.args[0]==uni :
                        coefun=np.prod(tmp.args[len(tmp.args)-1:])
                    if tmp.args[0]==uniplusun :
                        coefunplusun=np.prod(tmp.args[len(tmp.args)-1:])
                    if tmp.args[0]==unimoinsun :
                        coefunmoinsun=np.prod(tmp.args[len(tmp.args)-1:])
    return [coefunmoinsun*np.prod(coefficients), coefun*np.prod(coefficients),
     coefunplusun*np.prod(coefficients)]

def creation_matrice(coefs,maillage, pastemps, condition_de_bord = "Neumann"):
    """
    La condition de bord = Neumann ne sert à rien pour l'instant mais sera utile plus tard lorsqu'on voudra rajouter des options.

    La fonction prend en entrée la liste des coefficients retournés par la fonction calcul_coeffs, le pas de temps et le maillage.
    Elle retourne la matrice correspondante aux coefficients fournis (matrice constante pour l'instant).

    Si on arrive à décomposer les coefficients pour différentes variables (par exemple pour u et c dans Keller-Segel), on pourra 
    surement calculer les deux matrices indépendemment grâce à cette fonction.
    """
    taille=len(maillage)
    matrix = np.zeros((taille,taille))
    matrix[0][0]=coefs[0].subs({e : maillage[0], c : maillage[1],
         d : maillage[2], Deltat : pastemps, a : calcul_milieu(0,1),
          b : calcul_milieu(1,2) })
    matrix[0][1]=-coefs[0].subs({e : maillage[0], c : maillage[1],
         d : maillage[2], Deltat : pastemps, a : calcul_milieu(0,1),
          b : calcul_milieu(1,2) })
    for i in range(1,taille-1):
        matrix[i][i-1]=-coefs[0].subs({e : maillage[i-1], c : maillage[i],
         d : maillage[i+1], Deltat : pastemps, a : calcul_milieu(i-1,i),
          b : calcul_milieu(i,i+1) })
        matrix[i][i]=coefs[1].subs({e : maillage[i-1], c : maillage[i],
         d : maillage[i+1], Deltat : pastemps, a : calcul_milieu(i-1,i),
          b : calcul_milieu(i,i+1) })
        matrix[i][i+1]=-coefs[2].subs({e : maillage[i-1], c : maillage[i],
         d : maillage[i+1], Deltat : pastemps, a : calcul_milieu(i-1,i),
          b : calcul_milieu(i,i+1) })
    matrix[taille-1][taille-1]=coefs[0].subs({e : maillage[0], c : maillage[1],
         d : maillage[2], Deltat : pastemps, a : calcul_milieu(0,1),
          b : calcul_milieu(1,2) })
    matrix[taille-1][taille-2]=-coefs[0].subs({e : maillage[0], 
        c : maillage[1],
         d : maillage[2], Deltat : pastemps, a : calcul_milieu(0,1),
          b : calcul_milieu(1,2) })
    I=np.eye(len(maillage))
    return I-matrix

def animation_diffusion(maillage, matrice, tinit, tfinal, pasdetemps, u_init):
    """
    Fonction prenant en entrée le maillage ainsi que la matrice calculée avec creation_matrice,
     le temps initial, final et le pas de temps ainsi que le vecteur de coordonées initiales.
    La fonction ne retourne rien mais affiche l'évolution de la solution au cours du temps.
    """
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.set( ylim=(0, 1))
    ax.plot(maillage,u_init)
    while tinit<tfinal :
        u_init=np.dot(u_init,matrice)
        ax.clear()
        ax.set( ylim=(0, 1))
        ax.plot(maillage,u_init,'r')
        plt.pause(0.1)
        tinit+=pasdetemps


eq = equation_diffusion(0.5)
sp.preview(eq, viewer='file', filename='diffusion.png')
sp.pprint(eq)
expr = flux_dimension1(eq)
sp.preview(expr, viewer='file', filename='output1.png')
expr = reconstitution_differences_finies(expr)
sp.preview(expr, viewer='file', filename='volumes_finis.png')
print("#######################################################")
euler_explicite = schema_explicite(expr)
sp.preview(euler_explicite, viewer='file', filename='explicite.png')
maillage = maillage1D_regulier(0,1,0.02)
vecteurdescoeffcients = calcul_coeffs(euler_explicite)
mat = creation_matrice(vecteurdescoeffcients, maillage, 0.01)
print(mat)
u_init=np.ones((len(maillage)))
for i in range(0,len(u_init)):
    u_init[i]=0.5+0.5*random()
print(np.dot(mat,u_init))
animation_diffusion(maillage, mat, 0, 10, 0.01, u_init)