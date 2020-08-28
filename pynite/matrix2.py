import sympy as sp
import numpy as np
from pynite.tools import *
import matplotlib.pyplot as plt
import copy
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


################################################
#Version 2 pour que cela fonctionne avec un schéma implicite
################################################


def equation_factorisation(equation):
    """
    Parameters
    ----------
    equation: 
            equation after Euler scheme function
    Returns
    -------
        list of coefficients a,b,c (u^n+1 = a*u^n_i-1 + b*u^n_i + c*u^n_i+1) 
    """
    left = equation[0]
    right = equation[1]
    tmp = recomposition_basic_expression(left)
    left = decomposition_basic_expression(tmp)
    tmp  = recomposition_basic_expression(right)
    right = decomposition_basic_expression(tmp)
    print(right)
    coef = []
    for i in range(0,len(left)):
        for j in range(0,len(left[i])):
            if type(type(left[i][j]))==sp.core.function.UndefinedFunction :
                if left[i][j].args[1]==tn :
                    left[i].append(-1)
                    right.append(copy.deepcopy(left[i]))
                    left[i]=[0]
    for i in range(0,len(right)):
        for j in range(0,len(right[i])):
            if right[i]==[0]:
                0==0
            elif type(type(right[i][j]))==sp.core.function.UndefinedFunction :
                if right[i][j].args[1]==tnplusun :
                    right[i].append(-1)
                    left.append(copy.deepcopy(right[i]))
                    right[i]=[0]

    #A partir d'ici tout est à droite, il ne reste plus qu'à factoriser.
    #Cependant comme on ne connais pas forcément à l'avance la notation des variables (u et c), on doit les retrouver.
    variables_left = []
    symbole_variable_left = []
    for i in range(0,len(left)):
        for j in range(0,len(left[i])):
            if type(type(left[i][j]))==sp.core.function.UndefinedFunction :
                if left[i][j].args[0]!=xiplusundemi and left[i][j].args[0]!=ximoinsundemi :
                    if left[i][j] not in variables_left :
                        variables_left.append(left[i][j])
                    if type(left[i][j]) not in symbole_variable_left:
                        symbole_variable_left.append(type(left[i][j]))
    variables_right = []
    symbole_variable_right = copy.deepcopy(symbole_variable_left)
    for i in range(0,len(right)):
        for j in range(0,len(right[i])):
            if type(type(right[i][j]))==sp.core.function.UndefinedFunction :
                if right[i][j].args[0]!=xiplusundemi and right[i][j].args[0]!=ximoinsundemi :
                    if right[i][j] not in variables_right :
                        variables_right.append(right[i][j])
                    if type(right[i][j]) not in symbole_variable_right:
                        symbole_variable_right.append(type(right[i][j]))
    #On factorise à gauche
    args_mul = []
    coef=[]
    taille = len(symbole_variable_left)
    args_variables_left = [[0 for i in range(0,3)] for i in range(0,taille)]
    for i in range(0,len(left)) :
        test=0
        for j in range(0,len(left[i])):
            if left[i][j] in variables_left and test == 0:
                if left[i][j].args[0] == ximoinsun:
                    index = symbole_variable_left.index(type(left[i][j]))
                    args_variables_left[index][0] = np.prod(left[i][0:j]+left[i][j+1:len(left[i])]) + args_variables_left[index][0]
                    test=1
                elif left[i][j].args[0] == xi:
                    index = symbole_variable_left.index(type(left[i][j]))
                    args_variables_left[index][1] = np.prod(left[i][0:j]+left[i][j+1:len(left[i])]) + args_variables_left[index][1]
                    test=1
                elif left[i][j].args[0] == xiplusun:
                    index = symbole_variable_left.index(type(left[i][j]))
                    args_variables_left[index][2] = np.prod(left[i][0:j]+left[i][j+1:len(left[i])]) + args_variables_left[index][2]
                    test=1
    #On factorise à droite
    args_mul = []
    coef=[]
    taille = len(symbole_variable_right)
    args_variables_right = [[0 for i in range(0,3)] for i in range(0,taille)]
    for i in range(0,len(right)) :
        test=0
        for j in range(0,len(right[i])):
            if right[i][j] in variables_right and test==0:
                if right[i][j].args[0] == ximoinsun:
                    index = symbole_variable_right.index(type(right[i][j]))
                    args_variables_right[index][0] = np.prod(right[i][0:j]+right[i][j+1:len(right[i])]) + args_variables_right[index][0]
                    test=1
                elif right[i][j].args[0] == xi:
                    index = symbole_variable_right.index(type(right[i][j]))
                    args_variables_right[index][1] = np.prod(right[i][0:j]+right[i][j+1:len(right[i])]) + args_variables_right[index][1]
                    test=1
                elif right[i][j].args[0] == xiplusun:
                    index = symbole_variable_right.index(type(right[i][j]))
                    args_variables_right[index][2] = np.prod(right[i][0:j]+right[i][j+1:len(right[i])]) + args_variables_right[index][2]
                    test=1
    return ((args_variables_left, symbole_variable_left),(args_variables_right, symbole_variable_right))


###################################################
def matrix_expression(expression,mesh):
    """
    Parameters
    ----------
    expression: 
        list of lenght 2 : [sympy expression, symbols of variables] 
        (ex : [u_n+1 = 4*u_n + 3*u_n-1, u])
    mesh : list
        mesh 1D

    Returns
    -------
        list of matrix for the matrix form of the equation 
        (ex : [4*I, 3*I])
    """
    equations = expression[0]
    variables = expression[1]
    taille = len(mesh)
    matrix_sympy=[]
    matrix=[[[0 for i in range(0,taille)] 
            for j in range(0,taille)] for k in variables]
    for i in range(0,len(matrix)):
        for ligne in range(1,len(matrix[i])-1):
                matrix[i][ligne][ligne-1]=equations[i][0]
                matrix[i][ligne][ligne]=equations[i][1]
                matrix[i][ligne][ligne+1]=equations[i][2]
        matrix[i][0][1]=equations[i][2]
        matrix[i][taille-1][taille-2]=equations[i][0]
        if equations[i][0] != 0 :
            matrix[i][0][0]=equations[i][1]+equations[i][0]
            matrix[i][taille-1][taille-1]=equations[i][1]+equations[i][2]
        else :
            matrix[i][0][0]=equations[i][1]
            matrix[i][taille-1][taille-1]=equations[i][1]
        matrix_sympy.append(sp.Matrix(matrix[i]))

    return matrix_sympy

def matrix_equation(equation,mesh):
    """
    Parameters
    ----------
    equation: 
        sympy equation
    mesh : list
        mesh 1D
    Returns
    -------
        list of lenght 2 : [[matrix_left,symbols_left],[matrix_right,symbols_right]]
        matrix left and right are obtained with matrix_expression
    """
    left = equation[0]
    right = equation[1]
    matrix_left = matrix_expression(left,mesh)
    matrix_right = matrix_expression(right,mesh)
    return [[matrix_left,left[1]],[matrix_right,right[1]]]

def matrix_system(sys, mesh):
    """
    Parameters
    ----------
    sys: list 
        system represented as a list composed of sympy equations
    mesh : list
        mesh 1D
    Returns
    -------
        list of matricial representation of each equations in the system
    """
    matrix_sys=[]
    for lines in sys :
        matrix_sys.append(matrix_equation(lines,mesh))
    return matrix_sys

def matrix_sub(matrix, mesh, time_step):
    """
    Parameters
    ----------
    matrix: sympy Matrix
        symbolic matrix (one or more terms are symbols or functions)
    mesh : list
        mesh 1D
    time_step : float
    Returns
    -------
        Matrix with x_i,x_i+1,x_i-1,xi+1/2,xi-1/2,dt replaced by their float value
    """
    n = len(mesh)
    if matrix.is_diagonal()==True :
        for i in range(1,n-1):
            matrix[i*n + i] = matrix[i*n+i].subs({xi : mesh[i], Deltat : time_step,
                     xiplusun : mesh[i+1], ximoinsun : mesh[i-1],
                     xiplusundemi : middle(mesh[i],mesh[i+1]),
                     ximoinsundemi : middle(mesh[i],mesh[i-1])})
            matrix[0]=matrix[n+1]
            matrix[n**2-1]=matrix[n+1]
    else :
        matrix[0]=matrix[0].subs({xi : mesh[0], Deltat : time_step,
                     xiplusun : mesh[1],
                     xiplusundemi : middle(mesh[0],mesh[1]),
                     ximoinsundemi : mesh[0]-middle(mesh[0],mesh[1])})
        matrix[1]=matrix[1].subs({xi : mesh[0], Deltat : time_step,
                     xiplusun : mesh[1],
                     xiplusundemi : middle(mesh[0],mesh[1]),
                     ximoinsundemi : mesh[0]-middle(mesh[0],mesh[1])})
        for i in range(1,n-1):
            matrix[i*n + i] = matrix[i*n + i].subs({xi : mesh[i], Deltat : time_step,
                     xiplusun : mesh[i+1], ximoinsun : mesh[i-1],
                     xiplusundemi : middle(mesh[i],mesh[i+1]),
                     ximoinsundemi : middle(mesh[i],mesh[i-1])})
            matrix[i*n + i-1] = matrix[i*n + i-1].subs({xi : mesh[i], Deltat : time_step,
                     xiplusun : mesh[i+1], ximoinsun : mesh[i-1],
                     xiplusundemi : middle(mesh[i],mesh[i+1]),
                     ximoinsundemi : middle(mesh[i],mesh[i-1])})
            matrix[i*n + i+1] = matrix[i*n + i+1].subs({xi : mesh[i], Deltat : time_step,
                     xiplusun : mesh[i+1], ximoinsun : mesh[i-1],
                     xiplusundemi : middle(mesh[i],mesh[i+1]),
                     ximoinsundemi : middle(mesh[i],mesh[i-1])})
        matrix[n**2-1]=matrix[n**2-1].subs({xi : mesh[n-1], Deltat : time_step,
                     ximoinsun : mesh[n-2],
                     ximoinsundemi : middle(mesh[n-2],mesh[n-1]),
                     xiplusundemi : mesh[n-1]+mesh[n-1]-middle(mesh[n-1],mesh[n-2])})
        matrix[n**2-2]=matrix[n**2-2].subs({xi : mesh[n-1], Deltat : time_step,
                     ximoinsun : mesh[n-2],
                     ximoinsundemi : middle(mesh[n-2],mesh[n-1]),
                     xiplusundemi : mesh[n-1]+mesh[n-1]-middle(mesh[n-1],mesh[n-2])})
    return matrix



def matrix_subs_vector(matrix, vectors, symbols, mesh):
    """
    Parameters
    ----------
    matrix: sympy Matrix
        symbolic matrix (one or more terms are symbols or functions)
    vectors: list
        a list of floats values that represent U^n
    symbols : sympy symbol
        the symbol that represent the vector to replace
    mesh : list
        mesh 1D
    Returns
    -------
        Matrix with the symbol replaced by its float value
    """
    n = len(vectors[0])
    for cpt in range(0,len(symbols)) :
        matrix[0] = matrix[0].subs({symbols[cpt](mesh[0]-middle(mesh[0],mesh[1]),tn) : vectors[cpt][0],symbols[cpt](middle(mesh[0],mesh[1]),tn) : middle(vectors[cpt][0],vectors[cpt][1]),symbols[cpt](mesh[0],tn) : vectors[cpt][0], symbols[cpt](mesh[0],tnplusun) : vectors[cpt][0]})
        matrix[1] = matrix[1].subs({symbols[cpt](middle(mesh[0],mesh[1]),tn) : middle(vectors[cpt][0],vectors[cpt][1]),symbols[cpt](mesh[0]-middle(mesh[0],mesh[1]),tn) : vectors[cpt][0],symbols[cpt](mesh[0],tn) : vectors[cpt][0], symbols[cpt](mesh[0],tnplusun) : vectors[cpt][0]})
        for i in range(1,n-1):
            matrix[i*n + i] = matrix[i*n + i].subs({symbols[cpt](middle(mesh[i],mesh[i+1]),tn) : middle(vectors[cpt][i],vectors[cpt][i+1]), symbols[cpt](middle(mesh[i],mesh[i-1]),tn) : middle(vectors[cpt][i],vectors[cpt][i-1]), symbols[cpt](mesh[i],tn) : vectors[cpt][i], symbols[cpt](mesh[i],tnplusun) : vectors[cpt][i]})
            matrix[i*n + i-1] = matrix[i*n + i-1].subs({symbols[cpt](middle(mesh[i],mesh[i+1]),tn) : middle(vectors[cpt][i],vectors[cpt][i+1]), symbols[cpt](middle(mesh[i],mesh[i-1]),tn) : middle(vectors[cpt][i],vectors[cpt][i-1]), symbols[cpt](mesh[i],tn) : vectors[cpt][i], symbols[cpt](mesh[i],tnplusun) : vectors[cpt][i]})
            matrix[i*n + i+1] = matrix[i*n + i+1].subs({symbols[cpt](middle(mesh[i],mesh[i+1]),tn) : middle(vectors[cpt][i],vectors[cpt][i+1]), symbols[cpt](middle(mesh[i],mesh[i-1]),tn) : middle(vectors[cpt][i],vectors[cpt][i-1]), symbols[cpt](mesh[i],tn) : vectors[cpt][i], symbols[cpt](mesh[i],tnplusun) : vectors[cpt][i]})
        matrix[n**2-1] = matrix[n**2-1].subs({symbols[cpt](mesh[n-1]+mesh[n-1]-middle(mesh[n-1],mesh[n-2]),tn) : vectors[cpt][n-1],symbols[cpt](middle(mesh[n-1],mesh[n-2]),tn) : middle(vectors[cpt][n-1],vectors[cpt][n-2]), symbols[cpt](mesh[n-1],tn) : vectors[cpt][n-1], symbols[cpt](mesh[n-1],tnplusun) : vectors[cpt][n-1]})
        matrix[n**2-2] = matrix[n**2-2].subs({symbols[cpt](mesh[n-1]+mesh[n-1]-middle(mesh[n-1],mesh[n-2]),tn) : vectors[cpt][n-1],symbols[cpt](middle(mesh[n-1],mesh[n-2]),tn) : middle(vectors[cpt][n-1],vectors[cpt][n-2]), symbols[cpt](mesh[n-1],tn) : vectors[cpt][n-1], symbols[cpt](mesh[n-1],tnplusun) : vectors[cpt][n-1]})
    return matrix

def matrix_subs_sys(sys_matrix, vectors,mesh, time_step):
    """
    Parameters
    ----------
    sys_matrix: list
        representation of the matricial form of the system
    vectors : useless
    mesh : list
        mesh 1D
    time_step : float
    Returns
    -------
        system in matricial form with x_i,x_i+1,x_i-1,xi+1/2,xi-1/2,dt replaced by their float value
    """
    for lines in sys_matrix:
        left=lines[0]
        right=lines[1]
        matrix_left=left[0]
        symbols_left=left[1]
        matrix_right=right[0]
        symbols_right=right[1]
        matrix_left[0]=matrix_sub(matrix_left[0],mesh,time_step)
        for i in range(0,len(symbols_right)):
            matrix_right[i]=matrix_sub(matrix_right[i],mesh,time_step)
    return sys_matrix

def matrix_subs_equation(sys_matrix, vectors, mesh, time_step):
    """
    Parameters
    ----------
    sys_matrix: list
        representation of the matricial form of a system with a single equation in
    vectors : useless
    mesh : list
        mesh 1D
    time_step : float
    Returns
    -------
        system in matricial form with x_i,x_i+1,x_i-1,xi+1/2,xi-1/2,dt replaced by their float value
    """
    left=sys_matrix[0]
    right=sys_matrix[1]
    matrix_left=left[0]
    symbols_left=left[1]
    matrix_right=right[0]
    symbols_right=right[1]
    matrix_left[0]=matrix_sub(matrix_left[0],mesh,time_step)
    for i in range(0,len(symbols_right)):
        matrix_right[i]=matrix_sub(matrix_right[i],mesh,time_step)
    return sys_matrix

def animation_reworked(system,vectors,mesh,start,stop,time_step):
    """
    Parameters
    ----------
    system: list
        representation of the matricial form of the system 
    vectors : list of list of floats
        initial values of the vectors 
    mesh : list
        mesh 1D
    start : float
        initial time value
    stop : float
        stop time value
    time_step : float
        time step
    Returns
    -------
        None
    Prints
    -------
        The evolution of the system 
    """
    system_copy=copy.deepcopy(system)
    if len(vectors)!=1 :
        sys = matrix_subs_sys(system_copy,vectors,
        mesh ,time_step)
        symbols=sys[0][1][1]
        mat_right=[]
        mat_inv=[]
        for lines in sys :
            mat_left=lines[0][0][0]
            if mat_left.is_diagonal()==True :
                mat_inv.append(mat_left.inv())
            else :
                mat_inv.append(mat_left.inv(method="LU"))
            mat_right.append(lines[1][0])

        for i in range(0,len(mat_right)) :
            for j in range(0,len(mat_right[i])) :
                if mat_right[i][j].is_symbolic()==False:
                    mat_right[i][j] = mat_inv[i]*mat_right[i][j]
    else :
        sys = matrix_subs_equation(system_copy,vectors,mesh,time_step)
        symbols=sys[0][1]
        test=0
        mat_left=system_copy[0][0][0]
        mat_right=system_copy[1][0]
        mat_inv=[]
        if mat_left.is_symbolic()==True :
            test=1
            mat_inv.append(mat_left)
        elif mat_left.is_diagonal()==True :
            mat_inv.append(mat_left.inv())
        else :
            mat_inv.append(mat_left.inv(method="LU"))
        for i in range(0,len(mat_right)) :
            if mat_right[i].is_symbolic()==False:
                mat_right[i] = mat_inv[i]*mat_right[i]
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.set( ylim=(0, 1))
    ax.plot(mesh,vectors[0])
    plt.pause(0.01)
    while start < stop :
        start += time_step
        system_copy=copy.deepcopy(mat_right)
        if len(vectors)!=1 :
            for ligne in range(0,len(system_copy)) :
                for mat in range(0,len(system_copy[ligne])) :
                    if system_copy[ligne][mat].is_symbolic() == True :
                        system_copy[ligne][mat] = mat_inv[ligne]*matrix_subs_vector(system_copy[ligne][mat], vectors, symbols ,mesh)
        else :
            for mat in range(0,len(system_copy)):
                if system_copy[mat].is_symbolic()==True :
                    if test==1 :
                        #Pas encore bon pour l'equation de diffusion en implicite avec D(u) = D*u(x,t)
                        mat_inv_copy = copy.deepcopy(mat_inv[0])
                        inv = matrix_subs_vector(mat_inv_copy, vectors, symbols ,mesh).inv()
                        system_copy[mat] = inv*matrix_subs_vector(system_copy[mat], vectors, symbols ,mesh)
                    else :
                        system_copy[mat] = mat_inv[0]*matrix_subs_vector(system_copy[mat], vectors, symbols ,mesh)

        if len(vectors)!=1 :
            for i in range(0,len(vectors)):
                #On "remonte" le sytème en commençant par la dernière ligne
                tmp2=np.zeros(len(mesh))
                for j in range(0,len(system_copy[len(system_copy)-1-i])):
                    tmp2=np.dot(vectors[j],np.array(system_copy[len(sys)-1-i][j]).astype(np.float64))+tmp2
                vectors[len(vectors)-i-1] = tmp2
        else :
            vectors[0]=np.dot(vectors[0],np.array(system_copy[0]))

        ax.clear()
        ax.set(ylim=(0, 1))
        ax.plot(mesh,vectors[0],'r')
        plt.pause(0.01)

