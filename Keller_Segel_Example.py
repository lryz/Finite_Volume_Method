import pynite as pn
import sympy as sp

Un = sp.Symbol('U__n')
Unplusun = sp.Symbol('U__n+1')
Cn = sp.Symbol('C__n')
eq = pn.keller_segel(0.1,0.2,0.3,0.4,0.5)
sp.preview(eq, viewer='file', filename='keller_segel.png')
test=[]
for i in eq :
    test.append(pn.decomposition(i))

test2 = pn.flux_sys_1D(eq)
keller_segel_integre= []
for i in range(0,len(test2)):
    keller_segel_integre.append(pn.recomposition(test2[i]))
sp.preview(keller_segel_integre, viewer='file', filename='keller_segel_integre.png')
test=[]
for expr in keller_segel_integre :
    test.append(pn.finite_difference(pn.decomposition(expr)))

schema_explicit=[]
for i in range(0, len(test)):
    schema_explicit.append(pn.euler_forward(pn.decomposition(test[i])))
    schema_explicit[i] = pn.recomposition(schema_explicit[i])
    schema_explicit[i] = pn.decomposition(schema_explicit[i])
sp.preview(schema_explicit[0], viewer='file', filename='schema_explicite_ligne1.png')
sp.preview(schema_explicit[1], viewer='file', filename='schema_explicite_ligne2.png')
explicite_factorise=[]
for i in schema_explicit :
    explicite_factorise.append(pn.equation_factorisation(i))
sp.preview(explicite_factorise[0], viewer='file', filename='explicite_factorise_ligne1.png')
sp.preview(explicite_factorise[1], viewer='file', filename='explicite_factorise_ligne2.png')
sp.preview(explicite_factorise, viewer='file', filename='explicite_factorise.png')
maillage = pn.regular_mesh_1D(0,1,0.1)

matrices_l1 = pn.matrix_equation(explicite_factorise[0],maillage,0.1)
sp.preview(matrices_l1, viewer='file', filename='matrices_ligne_1.png')
matrices_l2 = pn.matrix_equation(explicite_factorise[1],maillage,0.1)
sp.preview(matrices_l2, viewer='file', filename='matrices_ligne_2.png')
sp.preview(sp.Eq(Unplusun * sp.eye(len(maillage)),matrices_l1[0]*Un + matrices_l1[1]*Cn), viewer='file', filename='matrice.png')