'''
'''




import pynite as pn
import sympy as sp
import numpy as np

u = sp.Function('u')
x = sp.Symbol('x')
t = sp.Symbol('t')
### Initialisation
## Definition of the constants
# Space :
dim = 1   # dimension
L = 1   # length of the domain
dx = 0.1  # grid step

# Temporal
dt = 0.0001 # time step
T = 10    # End time

# Model
D = 0.5 # diffusion coefficient

## Creation of the equation 
def test_diffusion(D):
    #D(u) = D*u(x,t)
    eq1 = sp.Eq(u(x,t).diff(t), D*u(x,t)*u(x,t).diff(x).diff(x))
    return eq1
eq = test_diffusion(D)

flux = pn.flux_sys_1D(eq)

diff_integrated = pn.recomposition(flux)

diff_reconstituted = pn.finite_difference(pn.decomposition(diff_integrated))


###Jusque la tout marche
sp.pprint(pn.decomposition(diff_reconstituted))
schema_explicit = pn.euler_forward(pn.decomposition(diff_reconstituted))

explicite_factorise = pn.equation_factorisation(schema_explicit)


mesh = pn.regular_mesh_1D(0,L,dx)

system = pn.matrix_equation(explicite_factorise,mesh)

u0=np.random.rand(len(mesh))
pn.animation_reworked(system, [u0],mesh,0,T, dt)