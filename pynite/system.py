import sympy as sp


x = sp.Symbol('x')
t = sp.Symbol('t')
u = sp.Function('u')
c = sp.Function('c')
phi = sp.Function('\phi')

def diffusion_equation(D):
    """
    La fonction equation_diffusion permet de générer un objet 
    de classe <class 'sympy.core.relational.Equality'>
    qui représente l'équation de diffusion.

    Elle prend en entrée le coefficient de diffusion D.
    """
    eq1 = sp.Eq(u(x,t).diff(t), D*u(x,t).diff(x).diff(x))
    return eq1


def keller_segel(k1,k2,kc,beta,alpha):
    """
        TODO: Explanation of the arguments. 
    """
    eq1 = sp.Eq(u(x,t).diff(t), k1*u(x,t).diff(x).diff(x) - k2*sp.Derivative(u(x,t)-u(x,t)**2,x)*c(x,t).diff(x).diff(x))
    eq2 = sp.Eq(c(x,t).diff(t), kc*c(x,t).diff(x).diff(x) - beta*c(x,t) + alpha*u(x,t)) 
    return [eq1,eq2]
