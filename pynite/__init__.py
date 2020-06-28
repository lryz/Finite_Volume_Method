import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from random import *
from pynite.matrix import *
from pynite.system import *
from pynite.mesh import *
from pynite.method import *
from pynite.tools import *
from pynite.interface_reconstitution import *

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
tn = sp.Symbol('t_n')
tnplusun = sp.Symbol('t_n+1')
Deltat = sp.Symbol('\Delta t')
uni = sp.Symbol('u__n_i')
unplusuni = sp.Symbol('u_i__n+1')
uniplusun = sp.Symbol('u_i+1__n')
unplusuniplusun = sp.Symbol('u_i+1__n')
unimoinsun = sp.Symbol('u_i-1__n')
U = sp.Symbol('U')
Un = sp.Symbol('U__n')
Unplusun = sp.Symbol('U__n+1')
Cn = sp.Symbol('C__n')
Cnplusun = sp.Symbol('C__n+1')
