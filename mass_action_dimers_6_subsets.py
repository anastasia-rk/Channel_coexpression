import matplotlib.pyplot as plt

from setup import *
import sympy as sym


# define symbolic variables
a1, a2, c1, c2, d1, d2, f1, f2, k1, k2 = sym.symbols('a1 a2 c1 c2 d1 d2 f1 f2 k1 k2')
z1, z2 = sym.symbols('z1 z2')
# define symbolic expressions: parts of mRNA that take part in the translation during normal function
u1 = f1 * (1 - c1) * z1 # mRNA1 that is colocalising with mRNA2 and not bound to ribosomes
u2 = f2 * (1 - c2) * z2 # mRNA2 that is colocalising with mRNA1 and not bound to ribosomes
w1 = (1 - f1) * (1 - c1) * z1 # mRNA1 that is not bound to ribosomes
w2 = (1 - f2) * (1 - c2) * z2 # mRNA2 that is not bound to ribosomes
y1 = (1 - d1) * c1 * z1 # mRNA1 that is bound to ribosomes and translated on its own
y2 = (1 - d2) * c2 * z2 # mRNA2 that is bound to ribosomes and translated on its own
x1 = d1 * c1 * z1 # mRNA1 that is bound to ribosomes and cotranslating with mRNA2
x2 = d2 * c2 * z2 # mRNA2 that is bound to ribosomes and cotranslating with mRNA1
# define symbolic equations: reduction in
eq0 = x1 - x2 # mRNA1 and mRNA2 that are cotranslating
eq0 = u1 - u2 # mRNA1 and mRNA2 that colocalising are not bound to ribosomes
eq1 = k1 * z1 - a1 * w1 - a1 * a2 * u1 - a1 * a2 * x1 - a1 * y1 # mRNA1 changes after gene silencing
eq2 = k2 * z2 - a2 * w2 - a1 * a2 * u2 - a1 * a2 * x2 - a2 * y2 # mRNA2 changes after gene silencing
sol_c1 = sym.solve(eq0, z1)[0]# solve for x1
sol_c2 = sym.solve(eq0, z2)[0]# solve for x2
# find functional form of a2 as a function of a1
a2_as_fn_of_k1a1 = sym.solve(eq1, a2)[0]# solve for a1
a2_as_fn_of_k2a1 = sym.solve(eq2, a2)[0]# solve for a2
eq_a2_is_a2 = a2_as_fn_of_k1a1 - a2_as_fn_of_k2a1 # the last two solutions give different expressions for a2
# we get the expression for a1 as a function of k1 and k2 and all the subpopulation parameters at normal function
a1_as_fn_of_k1k2 = sym.solve(eq_a2_is_a2, a1)# solve for a1
print(str(len(a1_as_fn_of_k1k2)) + ' unique solutions were found a1 as a function of k1 and k2')
# substitute the solution for a1 into the solution for a2
a2_as_fn_of_k1k2 = []
for i in range(len(a1_as_fn_of_k1k2)):
    a2_as_fn_of_k1k2.append(a2_as_fn_of_k2a1.subs(a1, a1_as_fn_of_k1k2[i]))
if len(a1_as_fn_of_k1k2) > 1:
    if a2_as_fn_of_k1k2[0] == a2_as_fn_of_k1k2[1]:
        print('the solutions for k2 are identical')
####################################################################################################################
# assign values to c1, c2, d1, d2, f1, f2, k1, k2 - from Eichel et.al. for hERG and SCN5A
# # how many transripts are undergoing translation
c1_val = 0.17
c2_val = 0.12 # from Catherine given to me in their lab
# # how many cells are in a pair with the other transcript (out of translating ones)
d1_val = 0.48
d2_val = 0.70 # did we grab this from Jameson et.al. Figure 2f? - seems so
# # how many transcripts are in a pair with the other transcript  (out of non-translating ones)
f1_val = 0.16/(1-c1_val) # worked out from the number of non-translating dimers and non-translating total u/(1-c)z or u/(u-w)
f2_val = 0.15/(1-c2_val)
# # kappas - by how much total number of transcripts reduced after gene silencing
k1_eichel = 0.54
k2_eichel = 0.59
# # the membrane conductance levels after gene silencing - derived from the current data
cond_frac_herg_eichel = 0.441
cond_frac_scn5a_eichel = 0.53
gene_labels = ['hERG', 'SCN5A']
paperName = 'Eichel et.al.'
####################################################################################################################
# assign values to c1, c2, d1, d2, k1, k2 - from Jameson et.al. for hERG and CACNA1C
# #how many transripts are undergoing translation
# c1_val = 0.17
# # c2_val = c1_val/4 # this is our assumption
# c2_val = 0.19 # from Catherine given to me in their lab
# # how many cells are in a pair with the other transcript (out of translating ones)
# d1_val = 0.345
# d2_val = 0.72
# # how many transcripts are in a pair with the other transcript  (out of non-translating ones)
# f1_val = 0.1
# f2_val = 0.1
# # kappas - by how much total number of transcripts reduced after gene silencing
# k1_eichel = 0.52
# k2_eichel = 0.56
# # the membrane conductance levels after gene silencing - derived from the current data
# cond_frac_herg_eichel = 0.441
# cond_frac_scn5a_eichel = 0.37
# gene_labels = ['hERG', 'CACNA1C']
# paperName = 'Jameson et.al.'
####################################################################################################################
# # trying stuff for Lis - hERG and RYR2
# # how many transripts are undergoing translation
# c1_val = 0.17
# c2_val = c1_val # we have to make an assumption for this one
# # how many cells are in a pair with the other transcript (out of translating ones)
# d1_val = 0.06 # from Figure 1c in the paper
# d2_val = 0.05 # from Figure 2f in the paper
# # how many transcripts are in a pair with the other transcript  (out of non-translating ones)
# f1_val = 0.05 # we have to make an assumption about this
# f2_val = 0.05 # we have to make an assumption about this
# # kappas - by how much total number of transcripts reduced after gene silencing
# k1_eichel = 0.54
# k2_eichel = 0.98
# # the membrane conductance levels after gene silencing - derived from the current data
# cond_frac_herg_eichel = 0.441
# cond_frac_scn5a_eichel = 0.70 # this is arbitrary, there is no conductance associated with RYR2
# gene_labels = ['hERG', 'RYR2']
# paperName = 'Jameson et.al.'
# ####################################################################################################################
# # trying stuff for Lis - hERG and LAMP1
# # how many transripts are undergoing translation
# c1_val = 0.17
# c2_val = c1_val*4 # we have to make an assumption for this one
# # how many cells are in a pair with the other transcript (out of translating ones)
# d1_val = 0.01 # from Figure 1c in the paper
# d2_val = 0.001 # from Figure 2f in the paper
# # how many transcripts are in a pair with the other transcript  (out of non-translating ones)
# f1_val = 0.01 # we have to make an assumption about this
# f2_val = 0.01 # we have to make an assumption about this
# # kappas - by how much total number of transcripts reduced after gene silencing
# k1_eichel = 0.54
# k2_eichel = 0.98
# # the membrane conductance levels after gene silencing - derived from the current data
# cond_frac_herg_eichel = 0.441
# cond_frac_scn5a_eichel = 0.70 # this is arbitrary, there is no conductance associated with RYR2
# gene_labels = ['hERG', 'LAMP1']
# paperName = 'Jameson et.al.'
####################################################################################################################
# substitute values into the solution for a1 and a2
a1_as_fn_of_k1k2_hERG = []
a2_as_fn_of_k1k2_SCN5A = []
for i in range(len(a1_as_fn_of_k1k2)):
    a1_as_fn_of_k1k2_hERG.append(a1_as_fn_of_k1k2[i].subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val), (f1, f1_val), (f2, f2_val)]))
    a2_as_fn_of_k1k2_SCN5A.append(a2_as_fn_of_k1k2[i].subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val), (f1, f1_val), (f2, f2_val)]))
# check if the found solutions are unique
if len(a1_as_fn_of_k1k2_hERG) > 1:
    if a1_as_fn_of_k1k2_hERG[0] == a1_as_fn_of_k1k2_hERG[1]:
        print('the solutions for k1 are identical')
if len(a2_as_fn_of_k1k2_SCN5A) > 1:
    if a2_as_fn_of_k1k2_SCN5A[0] == a2_as_fn_of_k1k2_SCN5A[1]:
        print('the solutions for k2 are identical')
####################################################################################################################
# plot the resultant surfaces for k1 and k2 changing between 0.5 and 1
k1_val = np.linspace(0.3, 0.9, 100)
k2_val = np.linspace(0.3, 0.9, 100)
a1_as_fn_of_k1k2_hERG_val = []
a2_as_fn_of_k1k2_SCN5A_val = []
for i in range(len(a1_as_fn_of_k1k2_hERG)):
    a1_as_fn_of_k1k2_hERG_val.append(sym.lambdify((k1, k2), a1_as_fn_of_k1k2_hERG[i], 'numpy'))
    a2_as_fn_of_k1k2_SCN5A_val.append(sym.lambdify((k1, k2), a2_as_fn_of_k1k2_SCN5A[i], 'numpy'))
####################################################################################################################
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
l1_as_fn_of_k1k2 = l1_as_fn_of_a1a2.subs([(a1, a1_as_fn_of_k1k2[0]), (a2, a2_as_fn_of_k1k2[0])])
l2_as_fn_of_k1k2 = l2_as_fn_of_a1a2.subs([(a1, a1_as_fn_of_k1k2[0]), (a2, a2_as_fn_of_k1k2[0])])
# substitute numerical values that are known into the expressions
l1_as_fn_of_k1k2_hERG = l1_as_fn_of_k1k2.subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val), (f1, f1_val), (f2, f2_val)])
l2_as_fn_of_k1k2_SCN5A = l2_as_fn_of_k1k2.subs([(c1, c1_val), (c2, c2_val), (d1, d1_val), (d2, d2_val), (f1, f1_val), (f2, f2_val)])
# turn expressions into functions - we will need these for the transform in other files
l1_as_fn_of_k1k2_hERG_val = sym.lambdify((k1, k2), l1_as_fn_of_k1k2_hERG, 'numpy')
l2_as_fn_of_k1k2_SCN5A_val = sym.lambdify((k1, k2), l2_as_fn_of_k1k2_SCN5A, 'numpy')

if __name__ == '__main__':
    # plot stuff
    ####################################################################################################################
    # plot the surfaces of alphas as functions of kappas
    fig = plt.figure(figsize=(16, 8))
    k1_mesh, k2_mesh = np.meshgrid(k1_val, k2_val)
    a1_mesh = a1_as_fn_of_k1k2_hERG_val[0](k1_mesh, k2_mesh)
    a2_mesh = a2_as_fn_of_k1k2_SCN5A_val[0](k1_mesh, k2_mesh)
    # compute the points where k1=0.5 and k2=0.5
    a1_experiment = a1_as_fn_of_k1k2_hERG_val[0](k1_eichel, k2_eichel)
    a2_experiment = a2_as_fn_of_k1k2_SCN5A_val[0](k1_eichel, k2_eichel)
    # create a constraint surface as function of k1_mesh and k2_mesh at z=1
    constraint_mesh1 = np.ones(k1_mesh.shape)
    constraint_mesh2 = np.zeros(k1_mesh.shape)
    a_surfs = [a1_mesh, a2_mesh]
    a_points = [a1_experiment, a2_experiment]
    print('After gene silencing in Eichel et.al. the following changes in mRNA levels were observed:')
    print('Total ' + gene_labels[0] + ' mRNA count: ' + str(k1_eichel) + ' z1')
    print('Total ' + gene_labels[1] + ' mRNA count: ' + str(k2_eichel) + ' z2')
    print(gene_labels[0] + ' free mRNA: ' + str(a1_experiment) + ' w1')
    print(gene_labels[1] + ' free mRNA: ' + str(a2_experiment) + ' w2')
    print(gene_labels[0] + ' co-translating mRNA: ' + str(a1_experiment * a2_experiment) + ' x1')
    print(gene_labels[1] + ' co-translating mRNA: ' + str(a1_experiment * a2_experiment) + ' x2')
    print(gene_labels[0] + ' translating mRNA: ' + str(a1_experiment) + ' y1')
    print(gene_labels[1] + ' translating mRNA: ' + str(a2_experiment) + ' y2')
    for i in range(len(gene_labels)):
        # plot the surfaces
        ax = fig.add_subplot(1, 2, i + 1, projection='3d')
        ax.plot_surface(k1_mesh, k2_mesh, a_surfs[i], cmap='magma', alpha=0.7)
        ax.plot_surface(k1_mesh, k2_mesh, constraint_mesh1, color='k', shade=False, alpha=0.1,
                        label=r'$\alpha_{' + str(i + 1) + '} = 1$')
        ax.plot_surface(k1_mesh, k2_mesh, constraint_mesh2, color='k', shade=False, alpha=0.1,
                        label=r'$\alpha_{' + str(i + 1) + '} = 0$')
        ax.scatter(0.5, 0.58, a_points[i], color='orange', s=25, label=paperName)
        ax.set_xlabel(r'$\kappa_1$')
        ax.set_ylabel(r'$\kappa_2$')
        ax.set_zlabel(r'$\alpha_{' + str(i + 1) + '}$')
        ax.legend(loc='upper right', fontsize=10)
        ax.set_title(gene_labels[i] + ': solution 1')
    # create surfaces
    a1_mesh = a1_as_fn_of_k1k2_hERG_val[1](k1_mesh, k2_mesh)
    a2_mesh = a2_as_fn_of_k1k2_SCN5A_val[1](k1_mesh, k2_mesh)
    # compute the points where k1=0.5 and k2=0.5
    # a1_experiment = a1_as_fn_of_k1k2_hERG_val[1](k1_eichel, k2_eichel)
    # a2_experiment = a2_as_fn_of_k1k2_SCN5A_val[1](k1_eichel, k2_eichel)
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
    #     ax.plot_surface(k1_mesh, k2_mesh, a_surfs[i], cmap='magma', alpha=0.5)
    #     ax.plot_surface(k1_mesh, k2_mesh, constraint_mesh, color='k', shade=False, alpha=0.5)
    #     ax.scatter(0.5, 0.58, a_points[i], color='r', s=25, label='hERG1a silencing')
    #     ax.set_xlabel(r'$\sigma_1$')
    #     ax.set_ylabel(r'$\sigma_2$')
    #     ax.set_zlabel(r'$\alpha_{' + str(i+1) + '}$')
    #     ax.set_title(gene_labels[i] + ': solution 2')
    # plt.tight_layout(pad=0.3, w_pad=0.5, h_pad=0.5)
    plt.savefig('Figures/' + gene_labels[0] + '_' + gene_labels[1] + '_alphas_as_function_of_kappas_correct.png')
    ####################################################################################################################
    ## plot the surfaces that show how lambdas depend on sigmas
    l1_mesh = l1_as_fn_of_k1k2_hERG_val(k1_mesh, k2_mesh)
    l2_mesh = l2_as_fn_of_k1k2_SCN5A_val(k1_mesh, k2_mesh)
    l_surfs = [l1_mesh, l2_mesh]
    # compute the points where k1=0.5 and k2=0.5
    l1_experiment = l1_as_fn_of_k1k2_hERG_val(k1_eichel, k2_eichel)
    l2_experiment = l2_as_fn_of_k1k2_SCN5A_val(k1_eichel, k2_eichel)
    l_points = [l1_experiment, l2_experiment]
    fig = plt.figure(figsize=(16, 8))
    for i in range(len(gene_labels)):
        # plot the surfaces
        ax = fig.add_subplot(1, 2, i + 1, projection='3d')
        ax.plot_surface(k1_mesh, k2_mesh, l_surfs[i], cmap='magma', alpha=0.7)
        ax.scatter(k1_eichel, k2_eichel, l_points[i], color='orange', s=25, label=paperName)
        ax.set_xlabel(r'$\kappa_1$')
        ax.set_ylabel(r'$\kappa_2$')
        ax.set_zlabel(r'$\lambda_{' + str(i + 1) + '}$')
        ax.legend(loc='upper right', fontsize=10)
        ax.set_title(gene_labels[i])
    # plt.tight_layout(pad=0.3, w_pad=1.8, h_pad=1.8)
    plt.savefig('Figures/' + gene_labels[0] + '_' + gene_labels[1] + '_lambdas_as_function_of_kappas_correct.png')
    print(gene_labels[0] + ' protein generaton: ' + str(l1_experiment) + ' p1')
    print(gene_labels[1] + ' protein generation: ' + str(l2_experiment) + ' p2')
    ####################################################################################################################
    ## generate a scatter plot of the mRNA levels and the protein levels
    covars = pickle.load(open('Pickles/gtex_covariances.pkl', 'rb'))
    Sigma_bivariate = covars[gene_labels[1]]  # np.array([[0.25, 0.25*(0.5/0.7)], [0.25*(0.5/0.7), 0.25]])
    print('transcripts to vary: ', gene_labels)
    print('Sigmas: ', Sigma_bivariate)
    sampled_from_normal = np.random.multivariate_normal(mean=[0, 0], cov=Sigma_bivariate, size=5000)
    sampled_silencing = sampled_from_normal[(sampled_from_normal < 0).all(axis=1)]
    transcript_levels = np.exp(sampled_silencing)
    # compute spearman correlation coefficient
    rho_transc, pval_transcr = sp.stats.spearmanr(transcript_levels[:,0], transcript_levels[:,1])
    # compute levels of free and translating mRNA
    a1_sample = a1_as_fn_of_k1k2_hERG_val[0](transcript_levels[:,0], transcript_levels[:,1])
    a2_sample = a2_as_fn_of_k1k2_SCN5A_val[0](transcript_levels[:,0], transcript_levels[:,1])
    a1a2_sample = a1_sample * a2_sample
    l1_sample = l1_as_fn_of_k1k2_hERG_val(transcript_levels[:,0], transcript_levels[:,1])
    l2_sample = l2_as_fn_of_k1k2_SCN5A_val(transcript_levels[:,0], transcript_levels[:,1])
    rho_prot, pval_prot = sp.stats.spearmanr(l1_sample, l2_sample)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes = axes.ravel()
    axes[0].scatter(transcript_levels[:,0],transcript_levels[:,1], color='k', s=5, alpha=0.15)
    axes[0].scatter(1, 1, color='orange', s=45, marker='X', label=paperName + ' scrambled shRNA')
    axes[0].scatter(k1_eichel, k2_eichel, color='orange', s=45, label=paperName +' hERG1b shRNA')
    axes[0].text(0.05, 0.95, r'$\rho = $' + str(round(rho_transc, 3)),
                 fontsize=14, transform=axes[0].transAxes)
    # axes[0].text(0.05, 0.9, r'$p$-val' + str(round(pval_transcr, 4)),
    #              fontsize=14, transform=axes[0].transAxes)
    axes[0].set_xlabel(r'$\kappa_1$')
    axes[0].set_ylabel(r'$\kappa_2$')
    axes[0].legend(loc='lower right', fontsize=10)
    axes[0].set_title('Total mRNA fraction')
    axes[1].scatter(a1_sample, a1a2_sample, color='k', s=5, alpha=0.15)
    axes[1].scatter(a1_experiment, a1_experiment*a2_experiment, color='magenta', s=45,  label='Transformed')
    axes[1].set_xlabel(r'$\alpha_1$ (free)')
    axes[1].set_ylabel(r'$\alpha_1\alpha_2$ (dimer with '+ gene_labels[1]+')')
    axes[1].legend(loc='lower right', fontsize=10)
    axes[1].set_title('translanting hERG mRNA fraction')
    axes[2].scatter(l1_sample, l2_sample, color='k', s=5, alpha=0.15)
    axes[2].scatter(l1_experiment, l2_experiment, color='magenta', s=45, label='Transformed')
    axes[2].scatter(1, 1, color='orange', s=45, marker='X',label=paperName + ' scrambled shRNA')
    axes[2].text(0.05, 0.95, r'$\rho = $' + str(round(rho_prot, 3)),
                 fontsize=14, transform=axes[2].transAxes)
    # axes[2].text(0.05, 0.9, r'$p$-val' + str(round(pval_prot, 4)),
    #              fontsize=14, transform=axes[2].transAxes)
    axes[2].scatter(cond_frac_herg_eichel, cond_frac_scn5a_eichel, color='orange', s=45, label=paperName + ' hERG1b shRNA')
    axes[2].set_xlabel(r'$\lambda_1$')
    axes[2].set_ylabel('$\lambda_2$')
    axes[2].legend(loc='lower right', fontsize=10)
    axes[2].set_title('Protein fraction')
    for ax in axes:
        ax.set_xticks(np.linspace(0.1, 1, 10))
        ax.set_yticks(np.linspace(0.1, 1, 10))
    plt.tight_layout()
    figName = 'Figures/' + gene_labels[0] + '_' + gene_labels[1] + '_transforms_gene_silencing_correct.png'
    plt.savefig(figName, dpi=600)

    ####################################################################################################################
    ## plot for the entire cloud
    transcript_levels = np.exp(sampled_from_normal)
    # compute spearman correlation coefficient
    rho_transc, pval_transcr = sp.stats.spearmanr(transcript_levels[:, 0], transcript_levels[:, 1])
    # compute levels of free and translating mRNA
    a1_sample = a1_as_fn_of_k1k2_hERG_val[0](transcript_levels[:, 0], transcript_levels[:, 1])
    a2_sample = a2_as_fn_of_k1k2_SCN5A_val[0](transcript_levels[:, 0], transcript_levels[:, 1])
    a1a2_sample = a1_sample * a2_sample
    l1_sample = l1_as_fn_of_k1k2_hERG_val(transcript_levels[:, 0], transcript_levels[:, 1])
    l2_sample = l2_as_fn_of_k1k2_SCN5A_val(transcript_levels[:, 0], transcript_levels[:, 1])
    rho_prot, pval_prot = sp.stats.spearmanr(l1_sample, l2_sample)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes = axes.ravel()
    axes[0].scatter(transcript_levels[:,0],transcript_levels[:,1], color='k', s=2, alpha=0.15)
    axes[0].scatter(1, 1, color='orange', s=15, marker='x', label=paperName + ' scrambled shRNA')
    axes[0].scatter(k1_eichel, k2_eichel, color='orange', s=15, label=paperName +' hERG1b shRNA')
    axes[0].text(0.05, 0.95, r'$\rho = $' + str(round(rho_transc, 3)),
                 fontsize=14, transform=axes[0].transAxes)
    # axes[0].text(0.05, 0.9, r'$p$-val' + str(round(pval_transcr, 4)),
    #              fontsize=14, transform=axes[0].transAxes)
    axes[0].set_xlabel(r'$\kappa_1$')
    axes[0].set_ylabel(r'$\kappa_2$')
    axes[0].legend(loc='lower right', fontsize=10)
    axes[0].set_title('Total mRNA fraction')
    axes[1].scatter(a1_sample, a1a2_sample, color='k', s=2, alpha=0.15)
    axes[1].scatter(a1_experiment, a1_experiment*a2_experiment, color='magenta', s=15,  label='Transformed')
    axes[1].set_xlabel(r'$\alpha_1$ (monomer)')
    axes[1].set_ylabel(r'$\alpha_1\alpha_2$ (dimer with '+ gene_labels[1]+')')
    axes[1].legend(loc='lower right', fontsize=10)
    axes[1].set_title('translanting hERG mRNA fraction')
    axes[2].scatter(l1_sample, l2_sample, color='k', s=2, alpha=0.15)
    axes[2].scatter(l1_experiment, l2_experiment, color='magenta', s=15, label='Transformed')
    axes[2].scatter(1, 1, color='orange', s=15, marker='x',label=paperName + ' scrambled shRNA')
    axes[2].text(0.05, 0.95, r'$\rho = $' + str(round(rho_prot, 3)),
                 fontsize=14, transform=axes[2].transAxes)
    # axes[2].text(0.05, 0.9, r'$p$-val' + str(round(pval_prot, 4)),
    #              fontsize=14, transform=axes[2].transAxes)
    axes[2].scatter(cond_frac_herg_eichel, cond_frac_scn5a_eichel, color='orange', s=15, label=paperName + ' hERG1b shRNA')
    axes[2].set_xlabel(r'$\lambda_1$')
    axes[2].set_ylabel('$\lambda_2$')
    axes[2].legend(loc='lower right', fontsize=10)
    axes[2].set_title('Protein fraction')
    for ax in axes:
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_xticks(np.linspace(0.1, 1, 10))
        # ax.set_yticks(np.linspace(0.1, 1, 10))
    plt.tight_layout()
    figName = 'Figures/' + gene_labels[0] + '_' + gene_labels[1] + '_transforms_all_correct.png'
    plt.savefig(figName, dpi=600)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs = axs.ravel()
    axs[0].scatter(transcript_levels[:, 0], a1a2_sample, color='k', s=2, alpha=0.15)
    axs[0].plot([0, 8], [0, 8], color='grey', alpha=0.5, linestyle='-',label=r'$\kappa=\alpha_1 \alpha_2$')
    axs[0].axvline(x=k1_eichel, color='grey', alpha=0.5, linestyle='--')
    axs[0].text(0.14, 0.75, r'$\kappa_1 = $' + str(round(k1_eichel, 4)),
                fontsize=12, transform=axs[0].transAxes)
    axs[0].axhline(y=a1_experiment*a2_experiment, color='grey', alpha=0.5, linestyle='--')
    axs[0].text(0.7, 0.08, r'$\alpha_1 \alpha_2 = $' + str(round(a1_experiment*a2_experiment, 4)),
                fontsize=12, transform=axs[0].transAxes)
    axs[0].scatter(1, 1, color='orange', s=20, marker='x', label=paperName + ' scrambled shRNA')
    axs[0].scatter(k1_eichel, a1_experiment*a2_experiment, color='magenta', s=20, label=paperName + ' hERG1b shRNA')
    axs[0].set_xlabel('Total transctips of hERG')
    axs[0].set_xscale('symlog', base=2, linthresh=0.125)
    axs[0].set_yscale('symlog', base=2, linthresh=0.125)
    axs[0].set_xticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
    axs[0].set_yticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
    axs[0].set_xlim([0.125, 8])
    axs[0].set_ylim([0.125, 8])
    axs[0].legend(bbox_to_anchor=(0., 1.02, 2.2, .102), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)
    axs[0].set_ylabel('In a dimer with ' + gene_labels[1])
    axs[1].scatter(transcript_levels[:, 1], a1a2_sample, color='k', s=2, alpha=0.15)
    axs[1].plot([0, 6], [0, 6], color='grey',alpha=0.5, linestyle='-',label=r'$\kappa_2=\alpha_1 \alpha_2$')
    axs[1].set_xlabel('Total transctips of '+ gene_labels[1])
    axs[1].axvline(x=k2_eichel, color='grey',alpha=0.5, linestyle='--')
    axs[1].text(0.14, 0.75, r'$\kappa_2 = $' + str(round(k2_eichel, 4)),
                fontsize=12, transform=axs[1].transAxes)
    axs[1].axhline(y=a1_experiment*a2_experiment, color='grey',alpha=0.5, linestyle='--')
    axs[1].text(0.7, 0.08, r'$\alpha_1 \alpha_2 = $' + str(round(a1_experiment*a2_experiment, 4)),
                fontsize=12, transform=axs[1].transAxes)
    axs[1].scatter(1, 1, color='orange', s=20, marker='x', label=paperName + ' scrambled shRNA')
    axs[1].scatter(k2_eichel, a1_experiment * a2_experiment, color='magenta', s=20, label=paperName + ' hERG1b shRNA')
    axs[1].set_xscale('symlog', base=2, linthresh=0.125)
    axs[1].set_yscale('symlog', base=2, linthresh=0.125)
    axs[1].set_xticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
    axs[1].set_yticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
    axs[1].set_xlim([0.125, 8])
    axs[1].set_ylim([0.125, 8])
    # axs[1].legend(loc='upper right', fontsize=12)
    axs[1].set_ylabel('In a dimer with ' + gene_labels[0])
    plt.tight_layout()
    figName = 'Figures/' + gene_labels[0] + '_' + gene_labels[1] + '_dimer_vs_total.png'
    plt.savefig(figName, dpi=600)
    ####################################################################################################################
    ## plot the Q-Q plots for the l1_sample and l2_sample values and their logs in a figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    axes = axes.ravel()
    sp.stats.probplot(l1_sample, plot=axes[0])
    axes[0].set_title('Q-Q plot for ' + gene_labels[0] + ' protein levels', fontsize=14)
    sp.stats.probplot(l2_sample, plot=axes[1])
    axes[1].set_title('Q-Q plot for ' + gene_labels[1] + ' protein levels', fontsize=14)
    sp.stats.probplot(np.log(l1_sample), plot=axes[2])
    axes[2].set_title('Q-Q plot for log ' + gene_labels[0] + ' protein levels', fontsize=14)
    sp.stats.probplot(np.log(l2_sample), plot=axes[3])
    axes[3].set_title('Q-Q plot for log ' + gene_labels[1] + ' protein levels', fontsize=14)
    plt.tight_layout()
    figName = 'Figures/' + gene_labels[0] + '_' + gene_labels[1] + 'protein_levels_QQ_plots.png'
    plt.savefig(figName, dpi=600)
    ####################################################################################################################
    # obtain the covariance matrix for log protein levels
    log_l1_sample = np.log(l1_sample)
    log_l2_sample = np.log(l2_sample)
    log_l1_sample_mean = np.mean(log_l1_sample)
    log_l2_sample_mean = np.mean(log_l2_sample)
    Sigma_log_proteins = np.cov(log_l1_sample, log_l2_sample)
    print('Sigma log proteins: ', Sigma_log_proteins)
    print('Mean log proteins: ', log_l1_sample_mean, log_l2_sample_mean)
    # save the covariance to a piclke file with gene name
    with open('Pickles/' + gene_labels[0] + '_' + gene_labels[1] + '_cotranslational_covariance.pkl', 'wb') as f:
        pickle.dump(Sigma_log_proteins, f)

