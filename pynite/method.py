import sympy as sp
x = sp.Symbol('x')
t = sp.Symbol('t')
u = sp.Function('u')
c = sp.Function('c')
phi = sp.Function('\phi')
tn = sp.Symbol('t_n')
tnplusun = sp.Symbol('t_n+1')
Deltat = sp.Symbol('\Delta t')

def euler_forward(equation):
    """
    Euler forward method
    Parameters
    ----------
    equation : 
        decomposed equation 
        
    Returns
    -------
        decomposed equation after replacing du(x,t)/dt by ((ux,tn+1)-u(x,tn))/deltat and all others t by tn
    """
    left = equation[0]
    for i in range(0,len(left)):

        for j in range(0,len(left[i])):
            if type(left[i][j]) == sp.core.function.Derivative :
                if left[i][j].args[1][0] == t :
                    #replacing du(x,t)/dt by ((ux,tn+1)-u(x,tn))/deltat
                    left[i][j] = (left[i][j].args[0].subs({t : tnplusun}) - left[i][j].args[0].subs({t : tn}))/Deltat
    right=equation[1]
    for i in range(0,len(right)):
        for j in range(0,len(right[i])):
            #replacing t by tn
            right[i][j]= right[i][j].subs({t : tn})
    return [left,right]