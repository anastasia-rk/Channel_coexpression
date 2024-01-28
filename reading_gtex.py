import cmapPy
import matplotlib.pyplot as plt
import pandas as pd
from cmapPy.pandasGEXpress.parse_gct import parse
from cmapPy.pandasGEXpress.write_gct import write
from setup import *

# main
if __name__ == '__main__':
    filename = '../GTEx_data/gene_tpm_heart_left_ventricle.gct'
    # test = parse(filename) # this does not work
    gtex_df = pd.read_csv(filename, sep='\t', skiprows=2, index_col=0)
    # find the gene names
    gene_names = gtex_df['Description'].tolist()
    # find the string CACNA1C in the list of gene names
    genes_of_interest = ['KCNH2', 'KCNQ1', 'CACNA1C', 'SCN5A', 'KCND3']
    # find all strings that contain string 'RIN3' in the list of gene names
    rins_of_interest = [s for s in gene_names if 'score' in s]

    # find the indices of the genes of interest
    indices_of_interest = []
    for gene in genes_of_interest:
        indices_of_interest.append(gene_names.index(gene))
    # only keep the genes of interest in the dataframe
    gtex_df = gtex_df.iloc[indices_of_interest]
    # export RIN numbers from the attributes table
    attributesFilename = '../GTEx_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
    attributes_df = pd.read_csv(attributesFilename, sep='\t', index_col=0)
    # associate RINs to the correct samples
    attributes_df = attributes_df[attributes_df['SMTSD'] == 'Heart - Left Ventricle']
    # drop all rows for which value in SMRIN column is below 7 or NaN
    attributes_df = attributes_df[attributes_df['SMRIN'].notna()]
    attributes_df = attributes_df[attributes_df['SMRIN'] > 7]
    high_rin_ids = attributes_df.index.tolist()
    all_ids = gtex_df.columns.tolist()[2:]
    # find if any o the high_rin_ids are in all_ids
    high_rin_id_matches = [s for s in high_rin_ids if s in all_ids]
    # only keep columns with that are in high_rin_id_matches
    gtex_high_rin_df = gtex_df[high_rin_id_matches].astype(float)

    ####################################################################################################################
    # plot samples on diffeent scales
    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    for iAx,ax in enumerate(axs):
        # plot values starting from column with index 3
        ax.scatter(gtex_df.iloc[0,3:], gtex_df.iloc[iAx+1,3:],s=5,marker='.',alpha=0.5,label='Heart - Left Ventricle')
        ax.scatter(gtex_high_rin_df.iloc[0, 3:], gtex_high_rin_df.iloc[iAx + 1, 3:], s=7, marker='o',color='darkorange',label='RIN > 7.0')
        ax.set(xlabel='KCNH2, TMP', ylabel=genes_of_interest[iAx+1]+',TMP')
        if iAx == 0:
            ax.legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('correlation_GTEx.png')

    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    for iAx,ax in enumerate(axs):
        # plot values starting from column with index 3
        ax.scatter(np.log2(np.array(gtex_df.iloc[0,3:].astype('float')) +1), np.log2(np.array(gtex_df.iloc[iAx+1,3:].astype('float'))+1),s=5,marker='.',alpha=0.5,color='magenta',label='Heart - Left Ventricle')
        ax.scatter(np.log2(gtex_high_rin_df.iloc[0, 3:]+1), np.log2(gtex_high_rin_df.iloc[iAx + 1, 3:]+1), s=7, marker='o',color='darkorange',label='RIN > 7.0')
        ax.set(xlabel='KCNH2, log(TMP+1)', ylabel=genes_of_interest[iAx+1]+',log(TMP+1)')
        if iAx == 0:
            ax.legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('correlation_GTEx_log2.png')

    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    for iAx,ax in enumerate(axs):
        # plot values starting from column with index 3
        ax.scatter(np.log(np.array(gtex_df.iloc[0,3:].astype('float')) +1), np.log(np.array(gtex_df.iloc[iAx+1,3:].astype('float'))+1),s=5,marker='.',color='magenta',alpha=0.5,label='Heart - Left Ventricle')
        ax.scatter(np.log(gtex_high_rin_df.iloc[0, 3:]+1), np.log(gtex_high_rin_df.iloc[iAx + 1, 3:]+1), s=7, marker='o',color='darkorange',label='RIN > 7.0')
        ax.set(xlabel='KCNH2, log(TMP+1)', ylabel=genes_of_interest[iAx+1]+',log(TMP+1)')
        if iAx == 0:
            ax.legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('correlation_GTEx_log.png')

    # for all entries in gtex_high_rin_df, calculate make Q-Q plot
    fig, axs = plt.subplots(3, 5, figsize=(20,12))
    axs = axs.ravel()
    for iAx, ax in enumerate(axs[:5]):
        # plot values starting from column with index 3
        sp.stats.probplot(gtex_high_rin_df.iloc[iAx, 3:].astype(float), plot=ax)
        ax.set(xlabel='Theoretical Quantiles', ylabel='Ordered Values', title=genes_of_interest[iAx]+',TMP')
    for iAx, ax in enumerate(axs[5:10]):
        # plot values starting from column with index 3
        sp.stats.probplot(np.log(gtex_high_rin_df.iloc[iAx, 3:].astype(float)+1), plot=ax)
        ax.set(xlabel='Theoretical Quantiles', ylabel='Ordered Values', title=genes_of_interest[iAx]+',log(TMP+1)')
    for iAx, ax in enumerate(axs[10:]):
        # plot values starting from column with index 3
        sp.stats.probplot(np.log(np.log(gtex_high_rin_df.iloc[iAx, 3:].astype(float)+1)), plot=ax)
        ax.set(xlabel='Theoretical Quantiles', ylabel='Ordered Values', title=genes_of_interest[iAx]+',log(log(TMP+1))')
    fig.tight_layout(pad=0.3)
    plt.savefig('QQ_GTEx.png')

    ####################################################################################################################
    # add one and convert to log scale for all entrie in gtex_high_rin_df
    gtex_high_rin_df_log = np.log(gtex_high_rin_df + 1)
    # simple sample covariances
    covariances = dict.fromkeys(genes_of_interest[1:])
    means = dict.fromkeys(genes_of_interest[:])
    covariances_log = dict.fromkeys(genes_of_interest[1:])
    means_log = dict.fromkeys(genes_of_interest[:])
    means['KCNH2'] = np.mean(gtex_high_rin_df.iloc[0,:])
    means_log['KCNH2'] = np.mean(gtex_high_rin_df_log.iloc[0,:])
    for i in range(1, len(genes_of_interest)):
        covariances[genes_of_interest[i]] = np.cov(gtex_high_rin_df.iloc[0,:], gtex_high_rin_df.iloc[i,:])
        means[genes_of_interest[i]] = np.mean(gtex_high_rin_df.iloc[i,:])
        covariances_log[genes_of_interest[i]] = np.cov(gtex_high_rin_df_log.iloc[0,:], gtex_high_rin_df_log.iloc[i,:])
        means_log[genes_of_interest[i]] = np.mean(gtex_high_rin_df_log.iloc[i,:])
    print(means)
    print(covariances)
    print(means_log)
    print(covariances_log)
    # save covariances_log to a pickle file

    with open('Pickles/gtex_covariances.pkl', 'wb') as f:
        pickle.dump(covariances_log, f)

    ####################################################################################################################
    # for each entry in covariances_log, calculate the correlation coefficient
    correlations_log = dict.fromkeys(genes_of_interest[1:])
    for key in covariances_log.keys():
        correlations_log[key] = covariances_log[key][0,1] / np.sqrt(covariances_log[key][0,0] * covariances_log[key][1,1])
    print(correlations_log)

    ####################################################################################################################
    # for each entry in covariances_log, plot a countour of the bivariate normal distribution with mean_log and covariances_log
    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    for ikey, key in enumerate(covariances_log.keys()):
        mean_bivar = [means_log['KCNH2'], means_log[key]]
        cov_bivar  = covariances_log[key]
        # create a mesh grid with x and y values bounded by +- 3 standard deviations
        x, y = np.meshgrid(np.linspace(mean_bivar[0] - 3 * np.sqrt(cov_bivar[0,0]), mean_bivar[0] + 3 * np.sqrt(cov_bivar[0,0]), 100),
                            np.linspace(mean_bivar[1] - 3 * np.sqrt(cov_bivar[1,1]), mean_bivar[1] + 3 * np.sqrt(cov_bivar[1,1]), 100))
        # compute the bivariate normal pdf on the mesh grid
        z = np.exp(sp.stats.multivariate_normal.pdf(np.dstack((x, y)), mean=mean_bivar, cov=cov_bivar))
        # normalize the pdf
        z = z / np.max(z)
        # plot the contour
        cnt = axs[ikey].contour(np.exp(x), np.exp(y), z, cmap='plasma',levels=6)
        axs[ikey].clabel(cnt, cnt.levels, inline = True, fontsize = 10)
        # plot x=y line in dashed black
        axs[ikey].plot(np.exp(x[0,:]), np.exp(x[0,:]), '--k',alpha=0.5,label='x=y')
        axs[ikey].set_xlim(0, 100)
        axs[ikey].set_ylim(0, 60)
        axs[ikey].set_xlabel('KCNH2, TPM')
        axs[ikey].set_ylabel(key+', TPM')
        axs[ikey].legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('contour_bivariate_decimal.png')

    # for each entry in covariances_log, create a bivariate normal distribution with mean_log and covariances_log
    # and sample from it
    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    samples_bivar = dict.fromkeys(genes_of_interest[1:])
    for ikey, key in enumerate(covariances_log.keys()):
        mean_bivar = [means_log['KCNH2'], means_log[key]]
        cov_bivar  = covariances_log[key]
        x = np.random.multivariate_normal(mean_bivar, cov_bivar, 1000)
        samples_bivar[key] = x
        axs[ikey].scatter(x[:,0], x[:,1], s=7, marker='.', alpha=0.7,color='k',label='Sampled from normal')
        # find the best fit regression line and the correlation coefficient
        slope, intercept, r, p, se = sp.stats.linregress(x[:,0], x[:,1])
        axs[ikey].plot(x[:, 0], x[:, 0] * slope + intercept, lw='1', color='k',alpha=0.7,
                       label=r'$\rho$ synthetic= {:.4f}'.format(r))
        axs[ikey].scatter(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1),
                   np.log(np.array(gtex_df.iloc[ikey + 1, 3:].astype('float')) + 1), s=5, marker='.',color='magenta', alpha=0.5,
                   label='Heart - Left Ventricle')
        slope, intercept, r, p, se = sp.stats.linregress(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1), np.log(np.array(gtex_df.iloc[ikey + 1, 3:].astype('float')) + 1))
        axs[ikey].plot(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1), np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1) * slope + intercept, lw='1', color='magenta',
                label=r'$\rho$ GTEx= {:.4f}'.format(r))
        axs[ikey].scatter(np.log(gtex_high_rin_df.iloc[0, 3:] + 1), np.log(gtex_high_rin_df.iloc[ikey + 1, 3:] + 1), s=7,
                   marker='o',color='darkorange', label='RIN > 7.0')
        slope, intercept, r, p, se = sp.stats.linregress(np.log(gtex_high_rin_df.iloc[0, 3:] + 1), np.log(gtex_high_rin_df.iloc[ikey + 1, 3:] + 1))
        axs[ikey].plot(np.log(gtex_high_rin_df.iloc[0, 3:] + 1), np.log(gtex_high_rin_df.iloc[0, 3:] + 1) * slope + intercept, lw='1', color='darkorange',
                label=r'$\rho$ (RIN>7.0)= {:.4f}'.format(r))
        axs[ikey].set_xlabel('KCNH2, log(TMP+1)')
        axs[ikey].set_ylabel(key+',log(TMP+1)')
        axs[ikey].set_xscale('symlog', base=2)
        axs[ikey].set_yscale('symlog', base=2)
        axs[ikey].legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('sampled_from_bivariate_log.png')

    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    for ikey, key in enumerate(covariances_log.keys()):
        x = samples_bivar[key]
        axs[ikey].scatter(np.exp(x[:,0]), np.exp(x[:,1]), s=7, marker='.', alpha=0.7, color='k',label='Sampled from normal')
        # find spearman correlation of columns in x and plot the line  corresponding to that correlation
        slope, intercept, r, p, se = sp.stats.linregress(np.exp(x[:,0]), np.exp(x[:,1]))
        axs[ikey].plot(np.exp(x[:,0]), np.exp(x[:,0])*slope + intercept, lw='1', color='k',alpha=0.7, label=r'$\rho$ synthetic= {:.4f}'.format(r))
        axs[ikey].scatter(gtex_df.iloc[0, 3:], gtex_df.iloc[ikey + 1, 3:], s=5, marker='.',color='magenta', alpha=0.5,
                   label='Heart - Left Ventricle')
        slope, intercept, r, p, se = sp.stats.linregress(gtex_df.iloc[0, 3:].astype(float), gtex_df.iloc[ikey + 1, 3:].astype(float))
        axs[ikey].plot(gtex_df.iloc[0, 3:].astype(float), gtex_df.iloc[0, 3:].astype(float) * slope + intercept, lw='1', color='magenta',
                label=r'$\rho$ GTEx= {:.4f}'.format(r))
        axs[ikey].scatter(gtex_high_rin_df.iloc[0, 3:].astype(float), gtex_high_rin_df.iloc[ikey + 1, 3:].astype(float), s=7, marker='o',color='darkorange', label='RIN > 7.0')
        slope, intercept, r, p, se = sp.stats.linregress(gtex_high_rin_df.iloc[0, 3:].astype(float), gtex_high_rin_df.iloc[ikey + 1, 3:].astype(float))
        axs[ikey].plot(gtex_high_rin_df.iloc[0, 3:].astype(float), gtex_high_rin_df.iloc[0, 3:].astype(float) * slope + intercept, lw='1', color='darkorange',
                label=r'$\rho$ (RIN>7.0)= {:.4f}'.format(r))
        axs[ikey].set_xlabel('KCNH2, TMP')
        axs[ikey].set_ylabel(key+',TMP')
        axs[ikey].legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('sampled_from_bivariate_decimal.png')

    ####################################################################################################################
    # sample from zero mean multivariate normal distribution with covariances_log
    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    samples_bivar = dict.fromkeys(genes_of_interest[1:])
    for ikey, key in enumerate(covariances_log.keys()):
        mean_bivar = [0, 0]
        cov_bivar  = covariances_log[key]
        x = np.random.multivariate_normal(mean_bivar, cov_bivar, 1000)
        samples_bivar[key] = x
        axs[ikey].scatter(x[:,0], x[:,1], s=7, marker='.', alpha=0.7,color='k',label='Sampled from normal')
        # find the best fit regression line and the correlation coefficient
        slope, intercept, r, p, se = sp.stats.linregress(x[:,0], x[:,1])
        axs[ikey].plot(x[:, 0], x[:, 0] * slope + intercept, lw='1', color='k',alpha=0.7,
                       label=r'$\rho$ synthetic= {:.4f}'.format(r))
        # axs[ikey].scatter(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1) - means_log['KCNH2'],
        #            np.log(np.array(gtex_df.iloc[ikey + 1, 3:].astype('float')) + 1) - means_log[key], s=5, marker='.',color='magenta', alpha=0.5,
        #            label='Heart - Left Ventricle')
        # slope, intercept, r, p, se = sp.stats.linregress(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1) - means_log['KCNH2'], np.log(np.array(gtex_df.iloc[ikey + 1, 3:].astype('float')) + 1) - means_log[key])
        # axs[ikey].plot(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1)  - means_log['KCNH2'], (np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1) - means_log['KCNH2'])* slope + intercept, lw='1', color='magenta',
        #         label=r'$\rho$ GTEx= {:.4f}'.format(r))
        axs[ikey].scatter(np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'], np.log(gtex_high_rin_df.iloc[ikey + 1, 3:] + 1) - means_log[key], s=7,
                   marker='o',color='darkorange', label='RIN > 7.0')
        slope, intercept, r, p, se = sp.stats.linregress(np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'], np.log(gtex_high_rin_df.iloc[ikey + 1, 3:] + 1)- means_log[key])
        axs[ikey].plot(np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'], (np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'])* slope + intercept, lw='1', color='darkorange',
                label=r'$\rho$ (RIN>7.0)= {:.4f}'.format(r))
        axs[ikey].set_xlabel('KCNH2, log(TMP+1)')
        axs[ikey].set_ylabel(key+',log(TMP+1)')
        axs[ikey].legend()
        # # set x and y axes to symlog scales
        # axs[ikey].set_xscale('symlog', base=2)
        # axs[ikey].set_yscale('symlog', base=2)
    fig.tight_layout(pad=0.3)
    plt.savefig('sampled_from_bivariate_zeromean_log.png')

    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    for ikey, key in enumerate(covariances_log.keys()):
        x = samples_bivar[key]
        axs[ikey].scatter(np.exp(x[:,0]), np.exp(x[:,1]), s=7, marker='.', alpha=0.7, color='k',label='Sampled from normal')
        # find spearman correlation of columns in x and plot the line  corresponding to that correlation
        slope, intercept, r, p, se = sp.stats.linregress(np.exp(x[:,0]), np.exp(x[:,1]))
        axs[ikey].plot(np.exp(x[:,0]), np.exp(x[:,0])*slope + intercept, lw='1', color='k',alpha=0.7, label=r'$\rho$ synthetic= {:.4f}'.format(r))
        # axs[ikey].scatter(gtex_df.iloc[0, 3:] / means['KCNH2'], gtex_df.iloc[ikey + 1, 3:] / means[key], s=5, marker='.',color='magenta', alpha=0.5,
        #            label='Heart - Left Ventricle')
        # slope, intercept, r, p, se = sp.stats.linregress(gtex_df.iloc[0, 3:].astype(float)/means['KCNH2'], gtex_df.iloc[ikey + 1, 3:].astype(float)/means[key])
        # axs[ikey].plot(gtex_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_df.iloc[0, 3:].astype(float) * slope / means['KCNH2']+ intercept, lw='1', color='magenta',
        #         label=r'$\rho$ GTEx= {:.4f}'.format(r))
        axs[ikey].scatter(gtex_high_rin_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_high_rin_df.iloc[ikey + 1, 3:].astype(float) / means[key], s=7, marker='o',color='darkorange', label='RIN > 7.0')
        slope, intercept, r, p, se = sp.stats.linregress(gtex_high_rin_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_high_rin_df.iloc[ikey + 1, 3:].astype(float) / means[key])
        axs[ikey].plot(gtex_high_rin_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_high_rin_df.iloc[0, 3:].astype(float) * slope / means['KCNH2'] + intercept, lw='1', color='darkorange',
                label=r'$\rho$ (RIN>7.0)= {:.4f}'.format(r))
        axs[ikey].set_xlabel('KCNH2, TMP')
        axs[ikey].set_ylabel(key+',TMP')
        axs[ikey].legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('sampled_from_bivariate_zeromean_decimal.png')

    ####################################################################################################################
    # reduce the covariance matrix a bit to reduce the effect of heavy tailes in the GTEx data
    covariances_log_reduced = dict.fromkeys(genes_of_interest[1:])
    for key in covariances_log.keys():
        covariances_log_reduced[key] = covariances_log[key] * 0.7
    ####################################################################################################################
    # sample from zero mean multivariate normal distribution with covariances_log
    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    samples_bivar = dict.fromkeys(genes_of_interest[1:])
    for ikey, key in enumerate(covariances_log_reduced.keys()):
        mean_bivar = [0, 0]
        cov_bivar  = covariances_log_reduced[key]
        x = np.random.multivariate_normal(mean_bivar, cov_bivar, 1000)
        samples_bivar[key] = x
        axs[ikey].scatter(x[:,0], x[:,1], s=7, marker='.', alpha=0.7,color='k',label='Sampled from normal')
        # find the best fit regression line and the correlation coefficient
        slope, intercept, r, p, se = sp.stats.linregress(x[:,0], x[:,1])
        axs[ikey].plot(x[:, 0], x[:, 0] * slope + intercept, lw='1', color='k',alpha=0.7,
                       label=r'$\rho$ synthetic= {:.4f}'.format(r))
        # axs[ikey].scatter(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1) - means_log['KCNH2'],
        #            np.log(np.array(gtex_df.iloc[ikey + 1, 3:].astype('float')) + 1) - means_log[key], s=5, marker='.',color='magenta', alpha=0.5,
        #            label='Heart - Left Ventricle')
        # slope, intercept, r, p, se = sp.stats.linregress(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1) - means_log['KCNH2'], np.log(np.array(gtex_df.iloc[ikey + 1, 3:].astype('float')) + 1) - means_log[key])
        # axs[ikey].plot(np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1)  - means_log['KCNH2'], (np.log(np.array(gtex_df.iloc[0, 3:].astype('float')) + 1) - means_log['KCNH2'])* slope + intercept, lw='1', color='magenta',
        #         label=r'$\rho$ GTEx= {:.4f}'.format(r))
        axs[ikey].scatter(np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'], np.log(gtex_high_rin_df.iloc[ikey + 1, 3:] + 1) - means_log[key], s=7,
                   marker='o',color='darkorange', label='RIN > 7.0')
        slope, intercept, r, p, se = sp.stats.linregress(np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'], np.log(gtex_high_rin_df.iloc[ikey + 1, 3:] + 1)- means_log[key])
        axs[ikey].plot(np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'], (np.log(gtex_high_rin_df.iloc[0, 3:] + 1) - means_log['KCNH2'])* slope + intercept, lw='1', color='darkorange',
                label=r'$\rho$ (RIN>7.0)= {:.4f}'.format(r))
        axs[ikey].set_xlabel('KCNH2, log(TMP+1)')
        axs[ikey].set_ylabel(key+',log(TMP+1)')
        axs[ikey].legend()
        # axs[ikey].set_xscale('symlog', base=2)
        # axs[ikey].set_yscale('symlog', base=2)
    fig.tight_layout(pad=0.3)
    plt.savefig('sampled_from_bivariate_zeromean_sclaed_log.png')

    fig, axs = plt.subplots(2, 2, figsize=(10,10))
    axs = axs.ravel()
    for ikey, key in enumerate(covariances_log.keys()):
        x = samples_bivar[key]
        axs[ikey].scatter(np.exp(x[:,0]), np.exp(x[:,1]), s=7, marker='.', alpha=0.7, color='k',label='Sampled from normal')
        # find spearman correlation of columns in x and plot the line  corresponding to that correlation
        slope, intercept, r, p, se = sp.stats.linregress(np.exp(x[:,0]), np.exp(x[:,1]))
        axs[ikey].scatter(np.exp(x[:,0]), np.exp(x[:,0])*slope + intercept, s=2, marker='.', color='k',alpha=0.7, label=r'$\rho$ synthetic= {:.4f}'.format(r))
        # axs[ikey].scatter(gtex_df.iloc[0, 3:] / means['KCNH2'], gtex_df.iloc[ikey + 1, 3:] / means[key], s=5, marker='.',color='magenta', alpha=0.5,
        #            label='Heart - Left Ventricle')
        # slope, intercept, r, p, se = sp.stats.linregress(gtex_df.iloc[0, 3:].astype(float)/means['KCNH2'], gtex_df.iloc[ikey + 1, 3:].astype(float)/means[key])
        # axs[ikey].plot(gtex_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_df.iloc[0, 3:].astype(float) * slope / means['KCNH2']+ intercept, lw='1', color='magenta',
        #         label=r'$\rho$ GTEx= {:.4f}'.format(r))
        axs[ikey].scatter(gtex_high_rin_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_high_rin_df.iloc[ikey + 1, 3:].astype(float) / means[key], s=7, marker='o',color='darkorange', label='RIN > 7.0')
        slope, intercept, r, p, se = sp.stats.linregress(gtex_high_rin_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_high_rin_df.iloc[ikey + 1, 3:].astype(float) / means[key])
        axs[ikey].plot(gtex_high_rin_df.iloc[0, 3:].astype(float) / means['KCNH2'], gtex_high_rin_df.iloc[0, 3:].astype(float) * slope / means['KCNH2'] + intercept, lw='1', color='darkorange',
                label=r'$\rho$ (RIN>7.0)= {:.4f}'.format(r))
        axs[ikey].set_xlabel('KCNH2, TMP')
        axs[ikey].set_ylabel(key+',TMP')
        axs[ikey].legend()
    fig.tight_layout(pad=0.3)
    plt.savefig('sampled_from_bivariate_zeromean_scaled_decimal.png')

    print('pause here')