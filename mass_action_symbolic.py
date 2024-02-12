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
    # # assign values to c1, c2, d1, d2, s1, s2 - from Eichel et.al. for hERG and SCN5A
    # c1_val = 0.17
    # c2_val = 0.18
    # d1_val = 0.48
    # d2_val = 0.70
    # s1_eichel = 0.5
    # s2_eichel = 0.58
    # cond_frac_herg_eichel = 0.411
    # cond_frac_scn5a_eichel = 0.619
    # gene_labels = ['hERG', 'SCN5A']
    # paperName = 'Eichel et.al.'
    ####################################################################################################################
    # assign values to c1, c2, d1, d2, s1, s2 - from Jameson et.al. for hERG and CACNA1C
    c1_val = 0.17
    c2_val = 0.17
    d1_val = 0.345
    d2_val = 0.72
    s1_eichel = 0.52
    s2_eichel = 0.56
    cond_frac_herg_eichel = 0.411
    cond_frac_scn5a_eichel = 0.32
    gene_labels = ['hERG', 'CACNA1C']
    paperName = 'Jameson et.al.'
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
    fig = plt.figure(figsize=(16,8))
    s1_mesh, s2_mesh = np.meshgrid(s1_val, s2_val)
    a1_mesh = a1_as_fn_of_s1s2_hERG_val[0](s1_mesh, s2_mesh)
    a2_mesh = a2_as_fn_of_s1s2_SCN5A_val[0](s1_mesh, s2_mesh)
    # compute the points where s1=0.5 and s2=0.5
    a1_experiment = a1_as_fn_of_s1s2_hERG_val[0](s1_eichel, s2_eichel)
    a2_experiment = a2_as_fn_of_s1s2_SCN5A_val[0](s1_eichel, s2_eichel)
    # create a constraint surface as function of s1_mesh and s2_mesh at z=1
    constraint_mesh1 = np.ones(s1_mesh.shape)
    constraint_mesh2 = np.zeros(s1_mesh.shape)
    a_surfs = [a1_mesh, a2_mesh]
    a_points = [a1_experiment, a2_experiment]
    print('After gene silencing in Eichel et.al. the following changes in mRNA levels were observed:')
    print('Total '+ gene_labels[0] +' mRNA count: ' + str(s1_eichel) + ' z1')
    print('Total '+ gene_labels[1] +' mRNA count: ' + str(s2_eichel) + ' z2')
    print(gene_labels[0] +' free mRNA: ' + str(a1_experiment) + ' w1')
    print(gene_labels[1] +' free mRNA: ' + str(a2_experiment) + ' w2')
    print(gene_labels[0] +' co-translating mRNA: ' + str(a1_experiment*a2_experiment) + ' x1')
    print(gene_labels[1] +' co-translating mRNA: ' + str(a1_experiment*a2_experiment) + ' x2')
    print(gene_labels[0] +' translating mRNA: ' + str(a1_experiment) + ' y1')
    print(gene_labels[1] +' translating mRNA: ' + str(a2_experiment) + ' y2')
    for i in range(len(gene_labels)):
        # plot the surfaces
        ax = fig.add_subplot(1, 2, i+1, projection='3d')
        ax.plot_surface(s1_mesh, s2_mesh, a_surfs[i], cmap='magma', alpha=0.7)
        ax.plot_surface(s1_mesh, s2_mesh, constraint_mesh1,color='k', shade=False, alpha=0.1, label=r'$\alpha_{' + str(i+1) + '} = 1$')
        ax.plot_surface(s1_mesh, s2_mesh, constraint_mesh2, color='k', shade=False, alpha=0.1, label=r'$\alpha_{' + str(i+1) + '} = 0$')
        ax.scatter(0.5, 0.58, a_points[i], color='orange', s=25,label=paperName)
        ax.set_xlabel(r'$\kappa_1$')
        ax.set_ylabel(r'$\kappa_2$')
        ax.set_zlabel(r'$\alpha_{' + str(i+1) + '}$')
        ax.legend(loc='upper right',fontsize=10)
        ax.set_title(gene_labels[i] + ': solution 1')
    # create surfaces
    a1_mesh = a1_as_fn_of_s1s2_hERG_val[1](s1_mesh, s2_mesh)
    a2_mesh = a2_as_fn_of_s1s2_SCN5A_val[1](s1_mesh, s2_mesh)
    # compute the points where s1=0.5 and s2=0.5
    # a1_experiment = a1_as_fn_of_s1s2_hERG_val[1](s1_eichel, s2_eichel)
    # a2_experiment = a2_as_fn_of_s1s2_SCN5A_val[1](s1_eichel, s2_eichel)
    # a_surfs = [a1_mesh, a2_mesh]
    # a_points = [a1_experiment, a2_experiment]
    # print('After gene silencing in Eichel et.al. the following changes in mRNA levels were observed:')
    # print('hERG1a free mRNA: ' + str(a1_experiment) + ' w1')
    # print('SCN5A free mRNA: ' + str(a2_experiment) + ' w2')
    # print('hERG1a cotranslating mRNA: ' + str(a1_experiment*a2_experiment) + ' x')
    # print('SCN5A cotranslating mRNA: ' + str(a1_experiment*a2_experiment) + ' x')
    # print('hERG1a translating mRNA: ' + str(a1_experiment) + ' y1')
    # print('SCN5A translating mRNA: ' + str(a2_experiment) + ' y2')
    # for i in range(len(gene_labels)):
    #     # plot the surfaces
    #     ax = fig.add_subplot(2, 2, i+3, projection='3d')
    #     ax.plot_surface(s1_mesh, s2_mesh, a_surfs[i], cmap='magma', alpha=0.5)
    #     ax.plot_surface(s1_mesh, s2_mesh, constraint_mesh, color='k', shade=False, alpha=0.5)
    #     ax.scatter(0.5, 0.58, a_points[i], color='r', s=25, label='hERG1a silencing')
    #     ax.set_xlabel(r'$\sigma_1$')
    #     ax.set_ylabel(r'$\sigma_2$')
    #     ax.set_zlabel(r'$\alpha_{' + str(i+1) + '}$')
    #     ax.set_title(gene_labels[i] + ': solution 2')
    # plt.tight_layout(pad=0.3, w_pad=0.5, h_pad=0.5)
    plt.savefig('Figures/'+ gene_labels[0]+'_'+gene_labels[1]+'_alphas_as_function_of_kappas.png')
    ####################################################################################################################
    # model of the transaltion process
    # define symbolic variables
    # is rates have no relationship
    # l1, l2, r_yp1, r_yp2, r_xp1, r_xp2, rd = sym.symbols('l1 l2 r_yp1 r_yp2 r_xp1 r_xp2 rd')
    ## if rates have a relationship
    l1, l2, r_yp1, r_yp2, rd = sym.symbols('l1 l2 r_yp1 r_yp2 rd')
    # assume that rates of translation are not affected by the association of the mRNAs
    r_xp1 = r_yp1
    r_xp2 = r_yp2
    # define symbolic expressions: sythesis of proteins from mRNA
    p1 = (r_xp1/rd)*x1 + (r_yp1/rd)*y1 # protein1 that is translated from mRNA1
    p2 = (r_xp2/rd)*x2 + (r_yp2/rd)*y2 # protein2 that is translated from mRNA2
    # establish the relationship between the fractions of generated proteins from the colocalisation results
    # # p1 = ((r_xp1/rd)*x1)/d1
    # p1 = ((r_yp1/rd)*y1)/(1-d1)
    # # p2 = ((r_xp2/rd)*x2)/d2
    # p2 = ((r_yp2/rd)*y2)/(1-d2)
    r_yp1 = (1 - d1) * p1 * rd / y1
    r_yp2 = (1 - d2) * p2 * rd / y2
    # define symbolic equations: reduction in protein1 and protein2 after gene silencing
    eq_translation1 = l1*p1 - (r_xp1/rd)*a1*a2*x1 - (r_yp1/rd)*a1*y1
    eq_translation2 = l2*p2 - (r_xp2/rd)*a1*a2*x2 - (r_yp2/rd)*a2*y2
    # solve for l1 and l2 to see how they depend on alphas and normal production rates
    l1_as_fn_of_a1a2 = sym.solve(eq_translation1, l1)[0]
    l2_as_fn_of_a1a2 = sym.solve(eq_translation2, l2)[0]
    # substitute expressions for a1 and a2 into the solutions for l1 and l2
    l1_as_fn_of_s1s2 = l1_as_fn_of_a1a2.subs([(a1, a1_as_fn_of_s1s2[0]), (a2, a2_as_fn_of_s1s2[0])])
    l2_as_fn_of_s1s2 = l2_as_fn_of_a1a2.subs([(a1, a1_as_fn_of_s1s2[0]), (a2, a2_as_fn_of_s1s2[0])])
    # substitute numerical values that are known into the expressions
    l1_as_fn_of_s1s2_hERG = l1_as_fn_of_s1s2.subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val)])
    l2_as_fn_of_s1s2_SCN5A = l2_as_fn_of_s1s2.subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val)])
    # turn expressions into functions
    l1_as_fn_of_s1s2_hERG_val = sym.lambdify((s1, s2), l1_as_fn_of_s1s2_hERG, 'numpy')
    l2_as_fn_of_s1s2_SCN5A_val = sym.lambdify((s1, s2), l2_as_fn_of_s1s2_SCN5A, 'numpy')
    ## plot the surfaces that show how lambdas depend on sigmas
    l1_mesh = l1_as_fn_of_s1s2_hERG_val(s1_mesh, s2_mesh)
    l2_mesh = l2_as_fn_of_s1s2_SCN5A_val(s1_mesh, s2_mesh)
    l_surfs = [l1_mesh, l2_mesh]
    # compute the points where s1=0.5 and s2=0.5
    l1_experiment = l1_as_fn_of_s1s2_hERG_val(s1_eichel, s2_eichel)
    l2_experiment = l2_as_fn_of_s1s2_SCN5A_val(s1_eichel, s2_eichel)
    l_points = [l1_experiment, l2_experiment]
    fig = plt.figure(figsize=(16, 8))
    for i in range(len(gene_labels)):
        # plot the surfaces
        ax = fig.add_subplot(1, 2, i+1, projection='3d')
        ax.plot_surface(s1_mesh, s2_mesh, l_surfs[i], cmap='magma', alpha=0.7)
        ax.scatter(s1_eichel, s2_eichel, l_points[i], color='orange', s=25,label=paperName)
        ax.set_xlabel(r'$\kappa_1$')
        ax.set_ylabel(r'$\kappa_2$')
        ax.set_zlabel(r'$\lambda_{' + str(i+1) + '}$')
        ax.legend(loc='upper right',fontsize=10)
        ax.set_title(gene_labels[i])
    plt.tight_layout(pad=0.3, w_pad=0.5, h_pad=0.5)
    plt.savefig('Figures/'+ gene_labels[0]+'_'+gene_labels[1]+'_lambdas_as_function_of_kappas.png')
    print(gene_labels[0] + ' protein generaton: ' + str(l1_experiment) + ' p1')
    print(gene_labels[1] + ' protein generation: ' + str(l2_experiment) + ' p2')
    print('pause here')
    ####################################################################################################################
    covars = pickle.load(open('Pickles/gtex_covariances.pkl', 'rb'))
    Sigma_bivariate = covars[gene_labels[1]]  # np.array([[0.25, 0.25*(0.5/0.7)], [0.25*(0.5/0.7), 0.25]])
    print('transcripts to vary: ', gene_labels)
    print('Sigmas: ', Sigma_bivariate)
    sampled_from_normal = np.random.multivariate_normal(mean=[0, 0], cov=Sigma_bivariate, size=1000)
    sampled_from_normal = sampled_from_normal[(sampled_from_normal <= 0).all(axis=1)]
    transcript_levels = np.exp(sampled_from_normal)
    # compute levels of free and translating mRNA
    a1_sample = a1_as_fn_of_s1s2_hERG_val[0](transcript_levels[:,0], transcript_levels[:,1])
    a2_sample = a2_as_fn_of_s1s2_SCN5A_val[0](transcript_levels[:,0], transcript_levels[:,1])
    a1a2_sample = a1_sample * a2_sample
    l1_sample = l1_as_fn_of_s1s2_hERG_val(transcript_levels[:,0], transcript_levels[:,1])
    l2_sample = l2_as_fn_of_s1s2_SCN5A_val(transcript_levels[:,0], transcript_levels[:,1])
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes = axes.ravel()
    axes[0].scatter(transcript_levels[:,0],transcript_levels[:,1], color='k', s=5, alpha=0.35)
    axes[0].scatter(s1_eichel, s2_eichel, color='orange', s=25, label=paperName)
    axes[0].set_xlabel(r'$\kappa_1$')
    axes[0].set_ylabel(r'$\kappa_2$')
    axes[0].legend(loc='lower right', fontsize=10)
    axes[0].set_title('Total mRNA fraction')
    axes[1].scatter(a1_sample, a1a2_sample, color='k', s=5, alpha=0.35)
    axes[1].scatter(a1_experiment, a1_experiment*a2_experiment, color='magenta', s=25,  label='Transformed')
    axes[1].set_xlabel(r'$\alpha_1$')
    axes[1].set_ylabel(r'$\alpha_1\alpha_2$')
    axes[1].legend(loc='lower right', fontsize=10)
    axes[1].set_title('Free mRNA fraction')
    axes[2].scatter(l1_sample, l2_sample, color='k', s=5, alpha=0.35)
    axes[2].scatter(l1_experiment, l2_experiment, color='magenta', s=25, label='Transformed')
    axes[2].scatter(cond_frac_herg_eichel, cond_frac_scn5a_eichel, color='orange', s=25, label=paperName)
    axes[2].set_xlabel(r'$\lambda_1$')
    axes[2].set_ylabel('$\lambda_2$')
    axes[2].legend(loc='lower right', fontsize=10)
    axes[2].set_title('Protein fraction')
    for ax in axes:
        ax.set_xticks(np.linspace(0.1, 1, 10))
        ax.set_yticks(np.linspace(0.1, 1, 10))
    plt.tight_layout()
    figName = 'Figures/' + gene_labels[0] + '_' + gene_labels[1] + '_transforms_gene_silencing.png'
    plt.savefig(figName, dpi=600)

