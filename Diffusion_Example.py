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
D = 0.5 # diffusion coefficient

## Creation of the equation 
eq = pn.diffusion_equation(D)
sp.pprint(eq)
flux = pn.flux_sys_1D(eq)

diff_integrated = pn.recomposition(flux)

diff_reconstituted = pn.finite_difference(pn.decomposition(diff_integrated))

###Jusque la tout marche
schema_implicit = pn.euler_backward(pn.decomposition(diff_reconstituted))

implicite_factorise = pn.equation_factorisation(schema_implicit)

mesh = pn.regular_mesh_1D(0,L,dx)

system = pn.matrix_equation(implicite_factorise,mesh)

u0=np.random.rand(len(mesh))
pn.animation_reworked(system, [u0],mesh,0,T, dt)