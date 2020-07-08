import sympy as sp


x = sp.Symbol('x')
t = sp.Symbol('t')
u = sp.Function('u')
c = sp.Function('c')
phi = sp.Function('\phi')

def diffusion_equation(D):
    """
    Parameters
    ----------
    D: float
           diffusion coefficient
        
    Returns
    -------
        diffusion equation
    """
    eq1 = sp.Eq(u(x,t).diff(t), D*u(x,t).diff(x).diff(x))
    return eq1


def keller_segel(k1,k2,kc,beta,alpha):
    """
    Parameters
    ----------
    k1: float
            diffusion coefficient
    k2: float
            coefficient
    kc: float
            coefficient
    beta: float
            coefficient
    alpha: float
            coefficient
        
    Returns
    -------
        diffusion equation
    """
    eq1 = sp.Eq(u(x,t).diff(t), k1*u(x,t).diff(x).diff(x) - k2*(u(x,t)*c(x,t).diff(x)-sp.Pow(u(x,t),2)*c(x,t).diff(x)).diff(x,evaluate=False))
    eq2 = sp.Eq(c(x,t).diff(t), kc*c(x,t).diff(x).diff(x) - beta*c(x,t) + alpha*u(x,t))  
    return [eq1,eq2]
