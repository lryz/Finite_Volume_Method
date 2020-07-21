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
    Cette fonction prend en entrée notre équation factorisée après un schéma (uniquement Euler explicite pour l'instant) 
    et retourne une liste contenant les coefficients (dans l'ordre) des u^n_i-1, u^n_i, u^n_i+1
    """
    left = equation[0]
    right = equation[1]
    tmp = recomposition_basic_expression(left)
    left = decomposition_basic_expression(tmp)
    tmp  = recomposition_basic_expression(right)
    right = decomposition_basic_expression(tmp)
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
            if type(type(right[i][j]))==sp.core.function.UndefinedFunction :
                if right[i][j].args[1]==tnplusun :
                    right[i].append(-1)
                    left.append(copy.deepcopy(right[i]))
                    right[i]=[0]

    #A partir d'ici tout est à droite, il ne reste plus qu'à factoriser.
    #Cependant comme on ne connais pas forcément à l'avance la notation des variables (u et c), on doit les retrouver.
    sp.preview([left,right], viewer='file', filename='tmp.png')
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
    args_mul = []
    coef=[]
    taille = len(symbole_variable_left)
    args_variables_left = [[0 for i in range(0,3)] for i in range(0,taille)]
    for i in range(0,len(left)) :
        for j in range(0,len(left[i])):
            if left[i][j] in variables_left :
                if left[i][j].args[0] == ximoinsun:
                    index = symbole_variable_left.index(type(left[i][j]))
                    args_variables_left[index][0] = np.prod(left[i][0:j]+left[i][j+1:len(left[i])]) + args_variables_left[index][0]
                elif left[i][j].args[0] == xi:
                    index = symbole_variable_left.index(type(left[i][j]))
                    args_variables_left[index][1] = np.prod(left[i][0:j]+left[i][j+1:len(left[i])]) + args_variables_left[index][1]
                elif left[i][j].args[0] == xiplusun:
                    index = symbole_variable_left.index(type(left[i][j]))
                    args_variables_left[index][2] = np.prod(left[i][0:j]+left[i][j+1:len(left[i])]) + args_variables_left[index][2]

    args_mul = []
    coef=[]
    taille = len(symbole_variable_right)
    args_variables_right = [[0 for i in range(0,3)] for i in range(0,taille)]
    for i in range(0,len(right)) :
        for j in range(0,len(right[i])):
            if right[i][j] in variables_right :
                if right[i][j].args[0] == ximoinsun:
                    index = symbole_variable_right.index(type(right[i][j]))
                    args_variables_right[index][0] = np.prod(right[i][0:j]+right[i][j+1:len(right[i])]) + args_variables_right[index][0]
                elif right[i][j].args[0] == xi:
                    index = symbole_variable_right.index(type(right[i][j]))
                    args_variables_right[index][1] = np.prod(right[i][0:j]+right[i][j+1:len(right[i])]) + args_variables_right[index][1]
                elif right[i][j].args[0] == xiplusun:
                    index = symbole_variable_right.index(type(right[i][j]))
                    args_variables_right[index][2] = np.prod(right[i][0:j]+right[i][j+1:len(right[i])]) + args_variables_right[index][2]
    return ((args_variables_left, symbole_variable_left),(args_variables_right, symbole_variable_right))


###################################################
def matrix_expression(expression,mesh):
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
            matrix[i][0][0]=-matrix[i][0][1]
            matrix[i][taille-1][taille-1]=-matrix[i][taille-1][taille-2]
        else :
            matrix[i][0][0]=equations[i][1]
            matrix[i][taille-1][taille-1]=equations[i][1]
        matrix_sympy.append(sp.Matrix(matrix[i]))

    return matrix_sympy

def matrix_equation(equation,mesh):
    left = equation[0]
    right = equation[1]
    matrix_left = matrix_expression(left,mesh)
    matrix_right = matrix_expression(right,mesh)
    return [[matrix_left,left[1]],[matrix_right,right[1]]]

def matrix_system(sys, mesh):
    matrix_sys=[]
    for lines in sys :
        matrix_sys.append(matrix_equation(lines,mesh))
    return matrix_sys

def matrix_sub(matrix, mesh, time_step):
    ######
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
                     ximoinsundemi : 0})
        matrix[1]=matrix[1].subs({xi : mesh[0], Deltat : time_step,
                     xiplusun : mesh[1],
                     xiplusundemi : middle(mesh[0],mesh[1]),
                     ximoinsundemi : 0})
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
        matrix[n**2-1]=matrix[n**2-1].subs({xi : mesh[n-2], Deltat : time_step,
                     ximoinsun : mesh[n-1],
                     ximoinsundemi : middle(mesh[n-2],mesh[n-1]),
                     xiplusundemi : 0})
        matrix[n**2-2]=matrix[n**2-1].subs({xi : mesh[n-2], Deltat : time_step,
                     ximoinsun : mesh[n-1],
                     ximoinsundemi : middle(mesh[n-2],mesh[n-1]),
                     xiplusundemi : 0})
    return matrix



def matrix_subs_vector(matrix, vectors, symbols, mesh):
    n = len(vectors[1])
    for cpt in range(0,len(symbols)) :
        matrix[0] = matrix[0].subs({symbols[cpt](middle(mesh[0],mesh[1]),tn) : middle(vectors[cpt][0],vectors[cpt][1])})
        matrix[1] = matrix[1].subs({symbols[cpt](middle(mesh[0],mesh[1]),tn) : middle(vectors[cpt][0],vectors[cpt][1])})
        for i in range(1,n-1):
            matrix[i*n + i] = matrix[i*n + i].subs({symbols[cpt](middle(mesh[i],mesh[i+1]),tn) : middle(vectors[cpt][i],vectors[cpt][i+1]), symbols[cpt](middle(mesh[i],mesh[i-1]),tn) : middle(vectors[cpt][i],vectors[cpt][i-1])})
            matrix[i*n + i-1] = matrix[i*n + i-1].subs({symbols[cpt](middle(mesh[i],mesh[i+1]),tn) : middle(vectors[cpt][i],vectors[cpt][i+1]), symbols[cpt](middle(mesh[i],mesh[i-1]),tn) : middle(vectors[cpt][i],vectors[cpt][i-1])})
            matrix[i*n + i+1] = matrix[i*n + i+1].subs({symbols[cpt](middle(mesh[i],mesh[i+1]),tn) : middle(vectors[cpt][i],vectors[cpt][i+1]), symbols[cpt](middle(mesh[i],mesh[i-1]),tn) : middle(vectors[cpt][i],vectors[cpt][i-1])})
        matrix[n**2-1] = matrix[n**2-1].subs({symbols[cpt](middle(mesh[n-1],mesh[n-2]),tn) : middle(vectors[cpt][n-1],vectors[cpt][n-2])})
        matrix[n**2-2] = matrix[n**2-2].subs({symbols[cpt](middle(mesh[n-2],mesh[n-1]),tn) : middle(vectors[cpt][n-1],vectors[cpt][n-2])})
    return matrix

def matrix_subs_sys(sys_matrix, vectors,mesh, time_step):
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


def animation_reworked(system,vectors,mesh,start,stop,time_step):
    system_copy=copy.deepcopy(system)
    sys = matrix_subs_sys(system_copy,vectors,
    mesh ,time_step)
    symbols=sys[0][1][1]
    print(symbols)
    mat_right=[]
    mat_inv=[]
    for lines in sys :
        mat_left=lines[0][0][0]
        mat_inv.append(mat_left.inv())
        mat_right.append(lines[1][0])
    for i in range(0,len(mat_right)) :
        for j in range(0,len(mat_right[i])) :
            if mat_right[i][j].is_symbolic()==False:
                mat_right[i][j] = mat_inv[i]*mat_right[i][j]
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.plot(mesh,vectors[0])
    plt.pause(0.1)
    ####Jusque la on est à peu près bon
    while start < stop :
        print(start)
        #print(mat_right[0][1])
        print(vectors)
        start += time_step
        system_copy=copy.deepcopy(mat_right)
        for ligne in range(0,len(system_copy)) :
            for mat in range(0,len(system_copy[ligne])) :
                if system_copy[ligne][mat].is_symbolic() == True :
                    system_copy[ligne][mat] = mat_inv[ligne]*matrix_subs_vector(system_copy[ligne][mat], vectors, symbols ,mesh)
        for i in range(0,len(vectors)):
            #On "remonte" le sytème en commençant par la dernière ligne
            tmp2=np.zeros(len(mesh))
            for j in range(0,len(system_copy[len(system_copy)-1-i])):
                tmp2=np.dot(vectors[j],np.array(system_copy[len(sys)-1-i][j]).astype(np.float64))+tmp2
            vectors[len(vectors)-i-1] = tmp2
        ax.clear()
        ax.plot(mesh,vectors[0],'r')
        plt.pause(1)

############################################
#J'essaye un nouveau truc
###########################################

def animation_reworked_bis(system_matrix,vectors,mesh,start,stop,time_step):
    system_copy=copy.deepcopy(system_matrix)
    sys = matrix_subs_sys(system_copy,vectors,mesh ,time_step)
    test=[]
    for i in range(0,len(sys)) :
        test.append([])
        for j in range(0,len(sys[i])) :
            if sys[i][j].is_symbolic() == True :
                test[i].append(True)
            else :
                test[i].append(False)
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.plot(mesh,vectors[0])
    ax.set( ylim=(-0.5, 2))
    plt.pause(0.1)
    while start < stop :
        print(start)
        start += time_step
        system_copy=copy.deepcopy(sys)
        for ligne in range(0,len(test)) :
            for mat in range(0,len(test[ligne])) :
                if test[ligne][mat] == True :
                    system_copy[ligne][mat] = matrix_subs_vector(system_copy[ligne][mat], vectors, symbols ,mesh)
        for i in range(0,len(vectors)):
            #On "remonte" le sytème en commençant par la dernière ligne
            tmp2=0
            for j in range(0,len(sys[len(sys)-1-i])):
                tmp2=np.dot(vectors[j],np.array(system_copy[len(sys)-1-i][j]).astype(np.float64))+tmp2
            vectors[len(vectors)-i-1] = tmp2
        ax.clear()
        ax.set( ylim=(-0.5, 2))
        ax.plot(mesh,vectors[0],'r')
        plt.pause(0.1)