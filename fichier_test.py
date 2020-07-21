import pynite as pn
import sympy as sp
import numpy as np
import copy 
import matplotlib.pyplot as plt
Un = sp.Symbol('U__n')
Unplusun = sp.Symbol('U__n+1')
ximoinsundemi = sp.Symbol('x_i-1/2')
xiplusundemi = sp.Symbol('x_i+1/2')
Cn = sp.Symbol('C__n')
u = sp.Function('u')
c = sp.Function('c')
x = sp.Symbol('x')
t = sp.Symbol('t')
tn = sp.Symbol('t_n')
def keller_segel(k1,k2,kc,beta,alpha):
    """
    """
    eq1 = sp.Eq(u(x,t).diff(t), k1*u(x,t).diff(x).diff(x) - k2*(u(x,t)*c(x,t).diff(x)-sp.Pow(u(x,t),2)*c(x,t).diff(x)).diff(x,evaluate=False))
    eq2 = sp.Eq(c(x,t).diff(t), kc*c(x,t).diff(x).diff(x) - beta*c(x,t) + alpha*u(x,t)) 
    return [eq1,eq2]

eq = keller_segel(1,1,1,1,1)

sp.preview(eq, viewer='file', filename='keller_segel.png')
pn.pre_subs_x(eq[0],xiplusundemi)
sp.pprint(eq)
test2 = pn.flux_sys_1D(eq)
keller_segel_integre= []
for i in range(0,len(test2)) :
    keller_segel_integre.append(pn.recomposition(test2[i]))
sp.preview(keller_segel_integre, viewer='file', filename='keller_segel_integre.png')
test=[]
for expr in keller_segel_integre :
    test.append(pn.finite_difference(pn.decomposition(expr)))
sp.preview(test, viewer='file', filename='differences_finies.png')
######################################
schema_explicit=[]
schema_implicit=[]
for i in range(0, len(test)):
    schema_explicit.append(pn.euler_forward(pn.decomposition(test[i])))
    schema_implicit.append(pn.euler_backward(pn.decomposition(test[i])))
    schema_explicit[i] = pn.recomposition(schema_explicit[i])
    schema_explicit[i] = pn.decomposition(schema_explicit[i])
sp.preview(schema_explicit[0], viewer='file', filename='schema_explicite_ligne1.png')
sp.preview(schema_explicit[1], viewer='file', filename='schema_explicite_ligne2.png')
sp.preview(schema_implicit[0], viewer='file', filename='schema_implicite_ligne1.png')
sp.preview(schema_implicit[1], viewer='file', filename='schema_implicite_ligne2.png')
factorise=[]
factorise.append(pn.equation_factorisation(schema_implicit[0]))
factorise.append(pn.equation_factorisation(schema_implicit[1]))

sp.preview(factorise[0], viewer='file', filename='facto_ligne1.png')
sp.preview(factorise[1], viewer='file', filename='facto_ligne2.png')

mesh = pn.regular_mesh_1D(0,1,0.1)
time_step=0.0001

system = pn.matrix_system(factorise,mesh)

u0=np.random.rand(len(mesh))*0.1 + 0.5
c0=np.random.rand(len(mesh))*0.1 + 0.5

pn.animation_reworked(system, [u0,c0],mesh,0,1, time_step)


