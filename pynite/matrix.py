import sympy as sp
import numpy as np
from pynite.tools import *
import matplotlib.pyplot as plt
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

def equation_factorisation(equation):
    """
    Cette fonction prend en entrée notre équation factorisée après un schéma (uniquement Euler explicite pour l'instant) 
    et retourne une liste contenant les coefficients (dans l'ordre) des u^n_i-1, u^n_i, u^n_i+1
    """
    left = equation[0]
    right = equation[1]
    coef = []
    for i in range(0,len(left)):
        for j in range(0,len(left[i])):
            if type(type(left[i][j]))==sp.core.function.UndefinedFunction :
                if left[i][j].args[1]==tn :
                    right.append([])
                    for k in range(0,len(left[i])):
                        right[len(right)-1].append(-left[i][k])
                    left[i]=[0]
                elif left[i][j].args[1]==tnplusun :
                    coef.append([])
                    for k in range(0,len(left[i])):
                        if k!=j:
                            coef[len(coef)-1].append(left[i][k])
                            left[i][k]=1
    sum_coef=0
    for i in range(0,len(coef)):
        sum_coef = sum_coef+np.prod(coef[i])
    for i in range(0,len(right)):
        right[i].append(sp.simplify(1/sum_coef))
    #A partir d'ici tout est à droite, il ne reste plus qu'à factoriser.
    #Cependant comme on ne connais pas forcément à l'avance la notation des variables (u et c), on doit les retrouver.
    variables = []
    symbole_variable = []
    for i in range(0,len(right)):
        for j in range(0,len(right[i])):
            if type(type(right[i][j]))==sp.core.function.UndefinedFunction :
                if right[i][j].args[0]!=xiplusundemi and right[i][j].args[0]!=ximoinsundemi :
                    if right[i][j] not in variables :
                        variables.append(right[i][j])
                    if type(right[i][j]) not in symbole_variable:
                        symbole_variable.append(type(right[i][j]))
    equation = recomposition([left,right])
    rhs_facto = equation.rhs.factor(variables)
    args_mul = []
    coef=[]
    taille = len(symbole_variable)
    args_variables = [[0 for i in range(0,3)] for i in range(0,taille)]
    test=[]
    decomposition_mul(rhs_facto,args_mul)
    for i in range(0,len(args_mul)):
        if type(args_mul[i])!=sp.core.add.Add :
            coef.append(args_mul[i])
        else :
            for k in range(0,len(args_mul[i].args)):
                test.append([])
                decomposition_mul(args_mul[i].args[k],test[k])
                for j in range(0,len(test[k])):
                    if type(test[k][j]) in symbole_variable :
                        index = symbole_variable.index(
                                type(test[k][j]))
                        if(test[k][j].args[0])==ximoinsun :
                            args_variables[index][0]=np.prod(coef)*test[k][0]
                        elif(test[k][j].args[0])==xi :
                            args_variables[index][1]=np.prod(coef)*test[k][0]
                        elif(test[k][j].args[0])==xiplusun :
                            args_variables[index][2]=np.prod(coef)*test[k][0]
    return (args_variables, symbole_variable)

def matrix_equation(explicite_factorise,maillage):
    equations = explicite_factorise[0]
    variables = explicite_factorise[1]
    taille = len(maillage)
    matrix_sympy=[]
    matrix=[[[0 for i in range(0,taille)] 
            for j in range(0,taille)] for k in variables]
    for i in range(0,len(matrix)):
        if equations[i][0] != 0 :
            matrix[i][0][0]=equations[i][0]
            matrix[i][taille-1][taille-1]=equations[i][0]
        else :
            matrix[i][0][0]=equations[i][1]
            matrix[i][taille-1][taille-1]=equations[i][1]
        matrix[i][0][1]=-equations[i][0]
        matrix[i][taille-1][taille-2]=-equations[i][0]
        for ligne in range(1,len(matrix[i])-1):
                matrix[i][ligne][ligne-1]=-equations[i][0]
                matrix[i][ligne][ligne]=equations[i][1]
                matrix[i][ligne][ligne+1]=-equations[i][2]
        if sp.Matrix(matrix[i]).is_diagonal():
            matrix_sympy.append(sp.Matrix(matrix[i]))
        else :
            matrix_sympy.append(sp.eye(taille) - sp.Matrix(matrix[i]))
    return matrix_sympy



####WORK IN PROGRESS

def animation_matrix(sys_matrix,vectors_init, mesh, start, end, time_step):
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.set( ylim=(0, 1))
    ax.plot(mesh,vectors_init[0])
    plt.pause(0.1)
    while start < end :
        print(start)
        start = start + time_step
        matrix_sys = calcul_matrices(sys_matrix,vectors_init, mesh, start, end, time_step)
        for i in range(0,len(vectors_init)):
            vectors_init[len(vectors_init)-1-i] = np.dot(matrix_sys[len(vectors_init)
            -1-i][0],vectors_init[0]) 
            + np.dot(matrix_sys[len(vectors_init)-1-i][1],vectors_init[1]) 
        ax.clear()

        ax.plot(mesh,vectors_init[0],'r')
        plt.pause(0.1)

def calcul_matrices(sys_matrix,vectors, mesh, start, end, time_step):
    ######
    #ça fonctionne pas encore jsp pk
    n = len(mesh)
    cpt=0
    for lines in sys_matrix :
        for matrix in lines :
            symbols=[]
            for tmp in range(0,len(matrix)):
                pre_symbols(matrix[tmp],symbols)
            for j in range(0,len(matrix[0 : n])):
                if len(symbols)>0 :
                    for tmp in range(0,len(symbols)):
                        matrix[j]=matrix[j].subs({symbols[tmp](xiplusundemi,tn)
                         : (vectors[tmp][1]+vectors[tmp][2])/2,symbols[tmp](ximoinsundemi,tn) : 
                         (vectors[tmp][0]+vectors[tmp][1])/2})
                matrix[j]=matrix[j].subs({xi : mesh[1], Deltat : time_step,
                     xiplusun : mesh[2], ximoinsun : mesh[0],
                     xiplusundemi : middle(mesh[1],mesh[2]),
                     ximoinsundemi : middle(mesh[1],mesh[0]) })
            for i in range(1,n-1):
                for j in range(0,len(matrix[i*n:(i+1)*n])) :
                    if len(symbols)>0 :
                        for tmp in range(0,len(symbols)):
                            matrix[i*n+j]=matrix[i*n+j].subs({symbols[tmp](
                                xiplusundemi,tn)
                         : (vectors[tmp][i]+vectors[tmp][i+1])/2,
                         symbols[tmp](ximoinsundemi,tn) : 
                         (vectors[tmp][i]+vectors[tmp][i-1])/2})
                    matrix[i*n+j]=matrix[i*n+j].subs({xi : mesh[i], Deltat : time_step,
                     xiplusun : mesh[i+1], ximoinsun : mesh[i-1],
                     xiplusundemi : middle(mesh[i],mesh[i+1]),
                     ximoinsundemi : middle(mesh[i],mesh[i-1])})
            for j in range(0,len(matrix[(n-1)*n : n*n])):
                if len(symbols)>0 :
                    for tmp in range(0,len(symbols)):
                        matrix[(n-1)*n + j]=matrix[(n-1)*n + j].subs(
                            {symbols[tmp](xiplusundemi,tn)
                         : (vectors[tmp][len(vectors[tmp])-1]+
                            vectors[tmp][len(vectors[tmp])-2])/2
                         ,symbols[tmp](ximoinsundemi,tn) : 0})
                         #(vectors[tmp][len(vectors[tmp])-3]+vectors[tmp][len(vectors[tmp])-2])/2})
                matrix[(n-1)*n + j]=matrix[(n-1)*n + j].subs({xi : mesh[n-2], Deltat : time_step,
                     ximoinsun : mesh[n-3],xiplusun : mesh[n-1],
                     xiplusundemi : middle(mesh[n-1],mesh[n-2]),
                     ximoinsundemi : middle(mesh[n-2],mesh[n-3]),})
            matrix = np.dot(matrix,sp.eye(n))
            cpt+=1
    return sys_matrix