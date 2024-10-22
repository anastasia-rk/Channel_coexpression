from setup import *
import matplotlib.lines as mlines
import corner
import operator

# main
if __name__ == '__main__':
    # read the single cell data for the rhythmonome genes
    filename = '../SingleCellGene_data/cardiomyocytes_rhythmonome_only.csv'
    single_cell_df = pd.read_csv(filename)
    # only read values for ventricular myocytes for this comparison
    single_cell_df_ventricular_only = single_cell_df[single_cell_df['Cluster'].str.contains('Ventricular')]
    # get unique cluster names from the remaining table
    cluster_names = single_cell_df['Cluster'].unique()
    # drop the raw for whic cluster contains '15' - this is a mix of right and and left ventricles
    single_cell_df_left_ventricle = single_cell_df_ventricular_only[~single_cell_df_ventricular_only['Cluster'].str.contains('15')]
    # get the gene names - all columns except the first two
    gene_names = single_cell_df.columns.tolist()[2:]
    # check if any genes from rhythmonome are missing in gene_names and print thse
    missing_genes = [gene for gene in genes_in_rhythmonome if gene not in gene_names]
    if missing_genes is not None:
        print(f'Missing genes in single cell data: {missing_genes}')
    ###########################################################################
    # for now let's focus on the ones we were modelling
    gene_names = ['KCNH2', 'KCNQ1', 'CACNA1C', 'SCN5A', 'KCND3', 'RYR2', 'LAMP1']
    single_cell_df_all_cardiomyocites = single_cell_df[gene_names].astype(float)
    single_cell_df_ventricular_only = single_cell_df_ventricular_only[gene_names].astype(float)
    single_cell_df_left_ventricle = single_cell_df_left_ventricle[gene_names].astype(float)
    # drop all raws that have all zero values
    single_cell_df_all_cardiomyocites = single_cell_df_all_cardiomyocites.loc[(single_cell_df_all_cardiomyocites != 0).any(axis=1)]
    single_cell_df_ventricular_only = single_cell_df_ventricular_only.loc[(single_cell_df_ventricular_only != 0).any(axis=1)]
    single_cell_df_left_ventricle = single_cell_df_left_ventricle.loc[(single_cell_df_left_ventricle != 0).any(axis=1)]
    # drop all raws that have at any columns that have the same values
    single_cell_df_all_clean = single_cell_df_all_cardiomyocites.copy()
    for iGene, gene1 in enumerate(gene_names[1:]):
        for jGene, gene2 in enumerate(gene_names[:iGene]):
            single_cell_df_all_clean = single_cell_df_all_clean[operator.ne(single_cell_df_all_clean[gene1], single_cell_df_all_clean[gene2])]
    single_cell_df_ventricular_only_clean = single_cell_df_ventricular_only.copy()
    for iGene, gene1 in enumerate(gene_names[1:]):
        for jGene, gene2 in enumerate(gene_names[:iGene]):
            single_cell_df_ventricular_only_clean = single_cell_df_ventricular_only_clean[operator.ne(single_cell_df_ventricular_only_clean[gene1], single_cell_df_ventricular_only_clean[gene2])]
    single_cell_df_left_ventricle_clean = single_cell_df_left_ventricle.copy()
    for iGene, gene1 in enumerate(gene_names[1:]):
        for jGene, gene2 in enumerate(gene_names[:iGene]):
            single_cell_df_left_ventricle_clean = single_cell_df_left_ventricle_clean[operator.ne(single_cell_df_left_ventricle_clean[gene1], single_cell_df_left_ventricle_clean[gene2])]
    # save the cleaned data into files
    single_cell_df_all_clean.to_csv('../SingleCellGene_data/6_genes_cardiomyocytes_clean.csv', index=False)
    single_cell_df_ventricular_only_clean.to_csv('../SingleCellGene_data/6_genes_ventricular_clean.csv', index=False)
    single_cell_df_left_ventricle_clean.to_csv('../SingleCellGene_data/6_genes_left_ventricle_clean.csv', index=False)
    ####################################################################################################################
    # read the GTEx data for the rhythmonome genes
    filename = '../GTEx_data/gtex_rhythmonome_genes.csv'
    gtex_df = pd.read_csv(filename)
    # read the GTEx data for the rhythmonome genes with high RIN
    filename = '../GTEx_data/gtex_rhythmonome_genes_high_RIN.csv'
    gtex_high_rin_df = pd.read_csv(filename)
    # transpose the GTEx dataframes
    gtex_df = gtex_df.transpose()
    gtex_high_rin_df = gtex_high_rin_df.transpose()
    # make the first row the header
    gtex_df.columns = gtex_df.iloc[0]
    gtex_high_rin_df.columns = gtex_high_rin_df.iloc[0]
    # drop the first two rows - these just contain names of genes
    gtex_df = gtex_df.drop(gtex_df.index[:2])
    gtex_high_rin_df = gtex_high_rin_df.drop(gtex_high_rin_df.index[:2])
    # order the columns of the GTEx dataframes to match the order of the single cell data
    gtex_df = gtex_df[gene_names].astype(float)
    gtex_high_rin_df = gtex_high_rin_df[gene_names].astype(float)
    ####################################################################################################################
    # for all genes but the first one, make scatterplots with the first one for each dataframes
    nGenes = len(gene_names)
    # make scatter plots
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.ravel()
    for iGene in range(nGenes-1):
        gene = gene_names[iGene+1]
        # single cell data
        axs[iGene].scatter(single_cell_df_all_cardiomyocites[gene_names[0]], single_cell_df_all_cardiomyocites[gene], s=1, color='grey',label='All mycoytes')
        axs[iGene].scatter(single_cell_df_ventricular_only[gene_names[0]], single_cell_df_ventricular_only[gene], s=1,color='purple',label='All ventricular')
        axs[iGene].scatter(single_cell_df_left_ventricle[gene_names[0]], single_cell_df_left_ventricle[gene], s=1,color='orange',label='Left ventricle')
        axs[iGene].legend(loc='lower right')
        axs[iGene].set_xlabel(gene_names[0] +', log(expr/total*10e4')
        axs[iGene].set_ylabel(gene + ', log(expr/total*10e4')
        axs[iGene].set_xlim([0, 6])
        axs[iGene].set_ylim([0, 6])
    fig.suptitle(f'n cells all: {single_cell_df_all_cardiomyocites.shape[0]}; n cells ventricular: {single_cell_df_ventricular_only.shape[0]}; n cells left ventricle: {single_cell_df_left_ventricle.shape[0]}')
    plt.tight_layout()
    plt.savefig(f'Figures/single_cell_scatter_all_vs_{gene_names[0]}_{gene}.png')
    plt.close()
    # make scatter plots - cleaned up data
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.ravel()
    for iGene in range(nGenes - 1):
        gene = gene_names[iGene + 1]
        # single cell data
        axs[iGene].scatter(single_cell_df_all_clean[gene_names[0]], single_cell_df_all_clean[gene],
                           s=1, color='grey',label='All mycoytes')
        axs[iGene].scatter(single_cell_df_ventricular_only_clean[gene_names[0]], single_cell_df_ventricular_only_clean[gene], s=1,
                           color='purple',label='All ventricular')
        axs[iGene].scatter(single_cell_df_left_ventricle_clean[gene_names[0]], single_cell_df_left_ventricle_clean[gene], s=1,
                           color='orange',label='Left ventricle')
        axs[iGene].legend(loc='lower right')
        axs[iGene].set_xlabel(gene_names[0] + ', log(expr/total*10e4)')
        axs[iGene].set_ylabel(gene + ', log(expr/total*10e4)')
        axs[iGene].set_xlim([0, 6])
        axs[iGene].set_ylim([0, 6])
    fig.suptitle(f'n cells all: {single_cell_df_all_clean.shape[0]}; n cells ventricular: {single_cell_df_ventricular_only_clean.shape[0]}; n cells left ventricle: {single_cell_df_left_ventricle_clean.shape[0]}')
    plt.tight_layout()
    plt.savefig(f'Figures/single_cell_scatter_all_vs_{gene_names[0]}_{gene}_clean.png')
    plt.close()
    # same plot for GTEx data
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.ravel()
    for iGene in range(nGenes-1):
        gene = gene_names[iGene+1]
        # single cell data
        axs[iGene].scatter(gtex_df[gene_names[0]], gtex_df[gene], s=1,label='All RIN')
        axs[iGene].scatter(gtex_high_rin_df[gene_names[0]], gtex_high_rin_df[gene], s=1,label='High RIN')
        axs[iGene].legend(loc='lower right')
        axs[iGene].set_xlabel(gene_names[0] +', TMP')
        axs[iGene].set_ylabel(gene +', TMP')
    plt.tight_layout()
    plt.savefig(f'Figures/GTEx_scatter_all_vs_{gene_names[0]}_{gene}.png')
    # repeat the plot of GTEx data with log of TPM +1 values
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.ravel()
    for iGene in range(nGenes-1):
        gene = gene_names[iGene+1]
        # single cell data
        axs[iGene].scatter(np.log(gtex_df[gene_names[0]]+1), np.log(gtex_df[gene]+1), s=1,label='All RIN')
        axs[iGene].scatter(np.log(gtex_high_rin_df[gene_names[0]]+1), np.log(gtex_high_rin_df[gene]+1), s=1,label='High RIN')
        axs[iGene].legend(loc='lower right')
        axs[iGene].set_xlabel(gene_names[0]+',log(TMP+1)')
        axs[iGene].set_ylabel(gene +', log(TMP+1)')
        axs[iGene].set_xlim([0, 6])
        axs[iGene].set_ylim([0, 6])
    plt.tight_layout()
    plt.savefig(f'Figures/GTEx_scatter_all_vs_{gene_names[0]}_{gene}_log.png')

    ####################################################################################################################
    # get summary statistics of log expression values for each gene
    # single cell data
    single_cell_summary = single_cell_df_all_cardiomyocites.describe()
    # GTEx data
    gtex_summary = gtex_df.describe()
    # GTEx data with high RIN
    gtex_high_rin_summary = gtex_high_rin_df.describe()
    print('pause here')



