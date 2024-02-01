from setup import *
import sympy as sym

if __name__ == '__main__':
    # define symbolic variables
    a1, a2, c1, c2, d1, d2, s1, s2 = sym.symbols('a1 a2 c1 c2 d1 d2 s1 s2')
    z1, z2 = sym.symbols('z1 z2')
    # define symbolic expressions: parts of mRNA that take part in the translation during normal function
    w1 = (1 - c1) * z1 # mRNA1 that is not bound to ribosomes
    w2 = (1 - c2) * z2 # mRNA2 that is not bound to ribosomes
    y1 = (1 - d1) * c1 * z1 # mRNA1 that is bound to ribosomes and translated on its own
    y2 = (1 - d2) * c2 * z2  # mRNA2 that is bound to ribosomes and translated on its own
    x1 = d1 * c1 * z1 # mRNA1 that is bound to ribosomes and cotranslating with mRNA2
    x2 = d2 * c2 * z2 # mRNA2 that is bound to ribosomes and cotranslating with mRNA1
    # define symbolic equations: reduction in
    eq0 = x1 - x2 # cotranslating mRNA1 and mRNA2 that are cotranslating
    eq1 = s1 * z1 - a1 * w1 - a1 * a2 * x1 - a1 * y1 # mRNA1 changes after gene silencing
    eq2 = s2 * z2 - a2 * w2 - a1 * a2 * x2 - a2 * y2 # mRNA2 changes after gene silencing
    sol_c1 = sym.solve(eq0, z1)[0]# solve for x1
    sol_c2 = sym.solve(eq0, z2)[0]# solve for x2
    # find functional form of a2 as a function of a1
    a2_as_fn_of_s1a1 = sym.solve(eq1, a2)[0]# solve for a1
    a2_as_fn_of_s2a1 = sym.solve(eq2, a2)[0]# solve for a2
    eq_a2_is_a2 = a2_as_fn_of_s1a1 - a2_as_fn_of_s2a1 # the last two solutions give different expressions for a2
    # we get the expression for a1 as a function of s1 and s2 and all the subpopulation parameters at normal function
    a1_as_fn_of_s1s2 = sym.solve(eq_a2_is_a2, a1)# solve for a1
    print(str(len(a1_as_fn_of_s1s2)) + ' unique solutions were found a1 as a function of s1 and s2')
    # substitute the solution for a1 into the solution for a2
    a2_as_fn_of_s1s2 = []
    for i in range(len(a1_as_fn_of_s1s2)):
        a2_as_fn_of_s1s2.append(a2_as_fn_of_s2a1.subs(a1, a1_as_fn_of_s1s2[i]))
    if len(a1_as_fn_of_s1s2) > 1:
        if a2_as_fn_of_s1s2[0] == a2_as_fn_of_s1s2[1]:
            print('the solutions for s2 are identical')
    ####################################################################################################################
    # assign values to c1, c2, d1, d2, s1, s2 - from Eichel et.al. for hERG and SCN5A
    c1_val = 0.17
    c2_val = 0.18
    d1_val = 0.48
    d2_val = 0.70
    s1_eichel = 0.5
    s2_eichel = 0.58
    # substitute values into the solution for a1 and a2
    a1_as_fn_of_s1s2_hERG = []
    a2_as_fn_of_s1s2_SCN5A = []
    for i in range(len(a1_as_fn_of_s1s2)):
        a1_as_fn_of_s1s2_hERG.append(a1_as_fn_of_s1s2[i].subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val)]))
        a2_as_fn_of_s1s2_SCN5A.append(a2_as_fn_of_s1s2[i].subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val)]))
    # check if the found solutions are unique
    if len(a1_as_fn_of_s1s2_hERG) > 1:
        if a1_as_fn_of_s1s2_hERG[0] == a1_as_fn_of_s1s2_hERG[1]:
            print('the solutions for s1 are identical')
    if len(a2_as_fn_of_s1s2_SCN5A) > 1:
        if a2_as_fn_of_s1s2_SCN5A[0] == a2_as_fn_of_s1s2_SCN5A[1]:
            print('the solutions for s2 are identical')
    ####################################################################################################################
    # plot the resultant surfaces for s1 and s2 changing between 0.5 and 1
    s1_val = np.linspace(0.3, 0.9, 100)
    s2_val = np.linspace(0.3, 0.9, 100)
    a1_as_fn_of_s1s2_hERG_val = []
    a2_as_fn_of_s1s2_SCN5A_val = []
    for i in range(len(a1_as_fn_of_s1s2_hERG)):
        a1_as_fn_of_s1s2_hERG_val.append(sym.lambdify((s1, s2), a1_as_fn_of_s1s2_hERG[i], 'numpy'))
        a2_as_fn_of_s1s2_SCN5A_val.append(sym.lambdify((s1, s2), a2_as_fn_of_s1s2_SCN5A[i], 'numpy'))
    ####################################################################################################################
    # plot the surfaces in two subplots
    fig = plt.figure(figsize=(16,16))
    gene_labels = ['hERG', 'SCN5A']
    s1_mesh, s2_mesh = np.meshgrid(s1_val, s2_val)
    a1_mesh = a1_as_fn_of_s1s2_hERG_val[0](s1_mesh, s2_mesh)
    a2_mesh = a2_as_fn_of_s1s2_SCN5A_val[0](s1_mesh, s2_mesh)
    # compute the points where s1=0.5 and s2=0.5
    a1_experiment = a1_as_fn_of_s1s2_hERG_val[0](s1_eichel, s2_eichel)
    a2_experiment = a2_as_fn_of_s1s2_SCN5A_val[0](s1_eichel, s2_eichel)
    # create a constraint surface as function of s1_mesh and s2_mesh at z=1
    constraint_mesh = np.ones(s1_mesh.shape)
    a_surfs = [a1_mesh, a2_mesh]
    a_points = [a1_experiment, a2_experiment]
    print('After gene silencing in Eichel et.al. the following changes in mRNA levels were observed:')
    print('Total hERG1a mRNA count: ' + str(s1_eichel) + ' z1')
    print('Total SCN5A mRNA count: ' + str(s2_eichel) + ' z2')
    print('hERG1a free mRNA: ' + str(a1_experiment) + ' w1')
    print('SCN5A free mRNA: ' + str(a2_experiment) + ' w2')
    print('hERG1a cotranslating mRNA: ' + str(a1_experiment*a2_experiment) + ' x1')
    print('SCN5A cotranslating mRNA: ' + str(a1_experiment*a2_experiment) + ' x2')
    print('hERG1a translating mRNA: ' + str(a1_experiment) + ' y1')
    print('SCN5A translating mRNA: ' + str(a2_experiment) + ' y2')
    for i in range(len(gene_labels)):
    # create a meshgrid of values for s1 and s2
        # plot the surfaces
        ax = fig.add_subplot(2, 2, i+1, projection='3d')
        ax.plot_surface(s1_mesh, s2_mesh, a_surfs[i], cmap='magma', alpha=0.5)
        ax.plot_surface(s1_mesh, s2_mesh, constraint_mesh,color='k', shade=False, alpha=0.5, label='Constraint')
        ax.scatter(0.5, 0.58, a_points[i], color='r', s=25,label='hERG1a silencing')
        ax.set_xlabel(r'$\sigma_1$')
        ax.set_ylabel(r'$\sigma_2$')
        ax.set_zlabel(r'$\alpha_{' + str(i+1) + '}$')
        ax.legend(loc='upper right',fontsize=10)
        ax.set_title(gene_labels[i] + ': solution 1')
    # create surfaces
    a1_mesh = a1_as_fn_of_s1s2_hERG_val[1](s1_mesh, s2_mesh)
    a2_mesh = a2_as_fn_of_s1s2_SCN5A_val[1](s1_mesh, s2_mesh)
    # compute the points where s1=0.5 and s2=0.5
    a1_experiment = a1_as_fn_of_s1s2_hERG_val[1](s1_eichel, s2_eichel)
    a2_experiment = a2_as_fn_of_s1s2_SCN5A_val[1](s1_eichel, s2_eichel)
    a_surfs = [a1_mesh, a2_mesh]
    a_points = [a1_experiment, a2_experiment]
    print('After gene silencing in Eichel et.al. the following changes in mRNA levels were observed:')
    print('hERG1a free mRNA: ' + str(a1_experiment) + ' w1')
    print('SCN5A free mRNA: ' + str(a2_experiment) + ' w2')
    print('hERG1a cotranslating mRNA: ' + str(a1_experiment*a2_experiment) + ' x')
    print('SCN5A cotranslating mRNA: ' + str(a1_experiment*a2_experiment) + ' x')
    print('hERG1a translating mRNA: ' + str(a1_experiment) + ' y1')
    print('SCN5A translating mRNA: ' + str(a2_experiment) + ' y2')
    for i in range(len(gene_labels)):
        # plot the surfaces
        ax = fig.add_subplot(2, 2, i+3, projection='3d')
        ax.plot_surface(s1_mesh, s2_mesh, a_surfs[i], cmap='magma', alpha=0.5)
        ax.plot_surface(s1_mesh, s2_mesh, constraint_mesh, color='k', shade=False, alpha=0.5)
        ax.scatter(0.5, 0.58, a_points[i], color='r', s=25, label='hERG1a silencing')
        ax.set_xlabel(r'$\sigma_1$')
        ax.set_ylabel(r'$\sigma_2$')
        ax.set_zlabel(r'$\alpha_{' + str(i+1) + '}$')
        ax.set_title(gene_labels[i] + ': solution 2')
    plt.tight_layout(pad=0.3, w_pad=0.5, h_pad=0.5)
    plt.savefig('alphas_as_function_of_sigmas.png')

    ####################################################################################################################
    # model of the transaltion process
    # define symbolic variables
    l1, l2, rt11, rt12, rt21, rt22, rd1, rd2 = sym.symbols('l1 l2 rt11 rt12 rt21 rt22 rd1 rd2')
    # define symbolic expressions: parts of mRNA that take part in the translation during normal function
    p1 = (rt11/rd1)*x1 + (rt12/rd1)*y1 # protein1 that is translated from mRNA1
    p2 = (rt21/rd2)*x2 + (rt22/rd2)*y2 # protein2 that is translated from mRNA2
    # define symbolic equations: reduction in protein1 and protein2 after gene silencing
    eq_translation1 = l1*p1 - (rt11/rd1)*a1*a2*x1 - (rt12/rd1)*a1*y1
    eq_translation2 = l2*p2 - (rt21/rd2)*a1*a2*x2 - (rt22/rd2)*a2*y2
    # solve for rt12 and rt22
    rt12_as_fn_of_rt11 = sym.solve(eq_translation1, rt12)[0]
    rt22_as_fn_of_rt21 = sym.solve(eq_translation2, rt22)[0]






    print('pause here')
