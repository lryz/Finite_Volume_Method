'''
'''




import pynite as pn
import sympy as sp
import numpy as np


### Initialisation
## Definition of the constants
# Space :
dim = 1   # dimension
L = 1   # length of the domain
dx = 0.02  # grid step

# Temporal
dt = 0.0001 # time step
T = 10    # End time

# Model
D = 1 # diffusion coefficient
Chi = 5 # Coeff of chemosensitivity

## Definition of unknowns at time n+1 and n
u = sp.Function('u')
c = sp.Function('c')
Un = sp.Symbol('U__n')
Unplusun = sp.Symbol('U__n+1')
Cn = sp.Symbol('C__n')

## Creation of the equation 
eq = pn.keller_segel(D,Chi,1,0.5,0.5)
sp.preview(eq, viewer='file', filename='keller_segel.png')
test=[]
for i in eq :
    test.append(pn.decomposition(i))

test2 = pn.flux_sys_1D(eq)
keller_segel_integre= []
for i in range(0,len(test2)):
    keller_segel_integre.append(pn.recomposition(test2[i]))
sp.preview(keller_segel_integre, viewer='file', filename='keller_segel_integre.png')
test=[]
for expr in keller_segel_integre :
    test.append(pn.finite_difference(pn.decomposition(expr)))

schema_implicit=[]
for i in range(0, len(test)):
    schema_implicit.append(pn.euler_backward(pn.decomposition(test[i])))
    schema_implicit[i] = pn.recomposition(schema_implicit[i])
    schema_implicit[i] = pn.decomposition(schema_implicit[i])
implicite_factorise=[]
for i in schema_implicit :
    implicite_factorise.append(pn.equation_factorisation(i))

mesh = pn.regular_mesh_1D(0,L,dx)

system = pn.matrix_system(implicite_factorise,mesh)
#sp.preview(matrices_l1, viewer='file', filename='mat1.png')
#sp.preview(matrices_l2, viewer='file', filename='mat2.png')
u0=np.random.rand(len(mesh))
c0=np.random.rand(len(mesh))

pn.animation_reworked(system, [u0,c0],mesh,0,T, dt)