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
    On suppose que l'équation est "triée" c'est à dire que la partie gauche contient les dérivées par rapport à t
    Cette fonction va renvoyer le schéma explicite sous forme d'une equation.
    On va ensuite pourvoir la décomposer
    Et on va donc ensuite pouvoir la passer sous forme matricielle
    """
    left = equation[0]
    for i in range(0,len(left)):

        for j in range(0,len(left[i])):
            if type(left[i][j]) == sp.core.function.Derivative :
                if left[i][j].args[1][0] == t :
                    left[i][j] = (left[i][j].args[0].subs({t : tnplusun}) - left[i][j].args[0].subs({t : tn}))/Deltat
    right=equation[1]
    for i in range(0,len(right)):
        for j in range(0,len(right[i])):
            right[i][j]= right[i][j].subs({t : tn})
    return [left,right]