import sympy as sp
import numpy as np

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

def middle(start,end):
    """returns the middle of an interval
    Parameters
    ----------
    start: float
        left bound of interval
    end: float
        right bound of interval

    Returns
    -------
    (start+end)/2:
        middle of interval
    """
    return (start+end)/2

def pre_xiplusundemi(expr,test_variable):
    """
    Parameters
    ----------
    expr: sympy.core expression
        any expression 
    test_variable: empty list
        
    Modifies
    -------
    True is added to test_variable if xi+1/2 is in the expression 
    """
    for arg in expr.args:
        if arg == xiplusundemi :
            test_variable.append(True)
        pre_xiplusundemi(arg,test_variable)

def pre_ximoinsundemi(expr,test_variable):
    """
    Parameters
    ----------
    expr: sympy.core expression
        any expression 
    test_variable: empty list
        
    Modifies
    -------
    test_variable : True is added to test_variable if xi-1/2 is 
    in the expression 
    """
    for arg in expr.args:
        if arg == ximoinsundemi :
            test_variable.append(True)
        pre_ximoinsundemi(arg,test_variable)

def pre_symbols(expr,symbols):
    """
    Parameters
    ----------
    expr: sympy.core expression
        any expression 
    symbols: empty list
        
    Modifies
    -------
    symbols : any function that is in the expression is added to symbols
    """
    for arg in expr.args:
        if type(type(arg)) == sp.core.function.UndefinedFunction :
            if type(arg) not in symbols :
                symbols.append(type(arg))
        pre_symbols(arg,symbols)

def decomposition_add(expr, ls) :
    """
    Parameters
    ----------
    expr: sympy.core.add.Add
        any addition
    ls: empty list
        
    Modifies
    -------
    ls : add to ls all terms of the addition 
    """
    if type(expr)==sp.core.add.Add :
        for arg in expr.args :
            ls.append(arg)
            decomposition_add(arg,ls)

def decomposition_mul(expr, ls):
    """
    Parameters
    ----------
    expr: sympy.core.mul.Mul
        any multiplication
    ls: empty list
        
    Modifies
    -------
    ls : add to ls all terms of the multiplication
    """
    if type(expr)==sp.core.mul.Mul :
        for arg in expr.args :
            ls.append(arg)
            decomposition_mul(arg, ls)
            
def decomposition(eq):
    """
    Parameters
    ----------
    eq: sympy.core 
        any equation
        
    Returns
    -------
    [mullistleft,mulistright]:
        a list of lists of lists that is interprated like this :
        [[a]+[[c]*[b]] = [[a]*[b]]+[c]]
    """
    eqleft = sp.expand(eq.lhs)
    eqright = sp.expand(eq.rhs)
    sumlistleft=[]
    sumlistright=[]
    mullistleft=[]
    mullistright=[]
    decomposition_add(eqleft,sumlistleft)
    if sumlistleft==[]:
        sumlistleft=[eqleft]
    decomposition_add(eqright,sumlistright)
    if sumlistright==[]:
        sumlistright=[eqright]
    for i in sumlistleft :
        mullistleft.append([])
    for i in sumlistright :
        mullistright.append([])
    for i in range(0,len(sumlistleft)):
        decomposition_mul(sumlistleft[i],mullistleft[i])
    for i in range(0,len(sumlistright)):
        decomposition_mul(sumlistright[i],mullistright[i])
    for i in range(0,len(mullistleft)):
        if mullistleft[i]==[]:
            mullistleft[i] = [sumlistleft[i]]
    for i in range(0,len(mullistright)):
        if mullistright[i]==[]:
            mullistright[i] = [sumlistright[i]]
    return [mullistleft,mullistright]

def recomposition(terms_list):
    """
    Parameters
    ----------
    terms_list: list 
        a list of lists of lists like this :
        [[a],[[c],[b]] = [[[a],[b]],[c]]
        
    Returns
    -------
    sympy.core.relational.Equality:
        reconstitues the equality based on terms_list :
        a+c*b = a*b+c
    """
    left = terms_list[0]
    right = terms_list[1]
    eqright=0
    for i in range(0,len(right)):
        cpt=0
        testleft=[]
        testright=[]
        constant=[]
        for j in range(0,len(right[i])):
            if type(right[i][j])==sp.core.add.Add :
                cpt=cpt+1
        if cpt>1 :
            for j in range(0,len(right[i])):
                if type(right[i][j])==sp.core.add.Add :
                    testleft.append(0)
                    for args in right[i][j].args :
                        test_variable=[]
                        pre_ximoinsundemi(args,test_variable)
                        if len(test_variable)>0 :
                            if test_variable[0]==True :
                                testleft[len(testleft)-1]=testleft[len(
                                    testleft)-1]+args
                if type(right[i][j])==sp.core.add.Add :
                    testright.append(0)
                    for args in right[i][j].args :
                        test_variable=[]
                        pre_xiplusundemi(args,test_variable)
                        if len(test_variable)>0 :
                            if test_variable[0]==True :
                                testright[len(testright)-1]=testright[len(
                                    testright)-1]+args
                if (type(right[i][j])==sp.core.numbers.Integer or
                 (type(right[i][j])==sp.core.numbers.Float)):
                    constant.append(right[i][j])
            eqright = eqright + np.prod(constant)*(np.prod(testright) - 
                np.prod(testleft))
        else :
            eqright = eqright + np.prod(right[i])
    eqleft=0
    for i in range(0,len(left)):
        cpt=0
        testleft=[]
        testright=[]
        constant=[]
        for j in range(0,len(left[i])):
            if type(left[i][j])==sp.core.add.Add :
                cpt=cpt+1
        if cpt>1 :
            for j in range(0,len(left[i])):
                if type(left[i][j])==sp.core.add.Add :
                    testleft.append(left[i][j].args[0])
                if type(left[i][j])==sp.core.add.Add :
                    testright.append(left[i][j].args[1])
                if (type(left[i][j])==sp.core.numbers.Integer or
                 (type(left[i][j])==sp.core.numbers.Float)):
                    constant.append(left[i][j])
            eqleft = eqleft + np.prod(constant)*(np.prod(testleft) -
             np.prod(testright))
        else :
            eqleft = eqleft + np.prod(left[i])
    return sp.Eq(sp.factor(eqleft),eqright)


def flux_eq_1D(terms_list):
    """
    La fonction permet d'intégrer sur x l'équation de diffusion de manière
     "Formelle".
    On va séparer les arguments de chaque côté de l'équation en prenant les
     listes des arguments.
    On est obligé de créer des listes bis car on ne peut pas assigner des 
    valeurs à des variables de type "Derivative".

    #WIP#
    Pour l'instant, on ne peut pas rajouter de terme d'addition dans 
    l'équation. On est obligé d'avoir des produits de chaque cotés.
    """
    eqleft = terms_list[0]
    eqright = terms_list[1]
    for i in range(0,len(eqright)):
        for j in range(len(eqright[i])):
            if  ((type(eqright[i][j])==sp.core.numbers.Float or 
                type(eqright[i][j])==sp.core.numbers.Integer)
                 and len(eqright[i])==1):
                #C'est le cas ou on a un réel en addition et pas en 
                #coefficient 
                #On doit donc l'intégrer sur la cellule Ki
                eqright[i][j]=eqright[i][j]*(xiplusundemi-ximoinsundemi)
            elif (type(eqright[i][j])== sp.core.function.Derivative):
                if eqright[i][j].args[1][0] == x :
                    if eqright[i][j].args[1][1] == 1 :
                        eqright[i][j] = (eqright[i][j].args[0].subs(
                            {x : xiplusundemi}) 
                        - eqright[i][j].args[0].subs({x : ximoinsundemi}))
                    else :
                        eqright[i][j] = (sp.Derivative(
                            eqright[i][j].args[0].subs({x : xiplusundemi})
                        ,(eqright[i][j].args[1])[0],(eqright[i][j].args[1])[1]
                        -1)
                        - sp.Derivative(eqright[i][j].args[0].subs(
                            {x : ximoinsundemi}),
                           (eqright[i][j].args[1])[0],
                           (eqright[i][j].args[1])[1]-1))
                elif eqright[i][j].args[1][0] == t :
                    eqright[i][j]=eqright[i][j].subs(
                        {x : xi})*(xiplusundemi-ximoinsundemi)
            elif type(type(eqright[i][j]==sp.core.function.UndefinedFunction)):
                eqright[i][j]=eqright[i][j].subs({x : xi})
    for i in range(0,len(eqleft)):
        for j in range(len(eqleft[i])):
            if  ((type(eqleft[i][j])==sp.core.numbers.Float or 
                type(eqleft[i][j])==sp.core.numbers.Integer)
                 and len(eqleft[i])==1):
                #C'est le cas ou on a un réel en addition et pas en 
                #coefficient 
                #On doit donc l'intégrer sur la cellule Ki
                eqleft[i][j]=eqleft[i][j]*(xiplusundemi-ximoinsundemi)
            elif (type(eqleft[i][j])== sp.core.function.Derivative):
                if eqleft[i][j].args[1][0] == x :
                    if eqleft[i][j].args[1][1] == 1 :
                        eqleft[i][j] = (eqleft[i][j].args[0].subs(
                            {x : xiplusundemi}) 
                        -eqleft[i][j].args[0].subs({x : ximoinsundemi}))
                    else :
                        eqleft[i][j] = (sp.Derivative(
                            eqleft[i][j].args[0].subs({x : xiplusundemi})
                        ,(eqleft[i][j].args[1])[0],(eqleft[i][j].args[1])[1]-1)
                        - sp.Derivative(eqleft[i][j].args[0].subs(
                            {x : ximoinsundemi}),
                           (eqleft[i][j].args[1])[0],
                           (eqleft[i][j].args[1])[1]-1))
                elif eqleft[i][j].args[1][0] == t :
                    eqleft[i][j]=eqleft[i][j].subs(
                        {x : xi})*(xiplusundemi-ximoinsundemi)
            elif type(type(eqleft[i][j]==sp.core.function.UndefinedFunction)):
                eqleft[i][j]=eqleft[i][j].subs({x : xi})
    return [eqleft,eqright]

def flux_sys_1D(sys):
    tmp=[]
    if type(sys)== list:
        test=[]
        for i in sys :
            test.append(decomposition(i))
        for i in test:
            tmp.append(flux_eq_1D(i))
    else :
        test = decomposition(sys)
        tmp = flux_dimension1_eq(test)
    return tmp

