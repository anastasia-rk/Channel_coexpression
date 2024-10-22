from setup import *
import matplotlib.lines as mlines
import corner

# main
if __name__ == '__main__':
    filename = '../GTEx_data/gene_tpm_heart_left_ventricle.gct'
    # test = parse(filename) # this does not work
    gtex_df = pd.read_csv(filename, sep='\t', skiprows=2, index_col=0)
    # find the gene names
    gene_names = gtex_df['Description'].tolist()
    # find the string CACNA1C in the list of gene names
    # genes_of_interest = ['KCNH2', 'KCNQ1', 'CACNA1C', 'SCN5A', 'KCND3']
    genes_of_interest = genes_in_rhythmonome + ['LAMP1']
    print(f'Number of genes of interest: {len(genes_of_interest)}')
    # currents_of_interest = ['IKr', 'IKs', 'ICaL', 'INa', 'Ito']
    # find all strings that contain string 'RIN3' in the list of gene names
    rins_of_interest = [s for s in gene_names if 'score' in s]

    # find the indices of the genes of interest
    indices_of_interest = []
    for gene in genes_of_interest:
        # chek if gene is in gene_names
        if gene in gene_names:
            indices_of_interest.append(gene_names.index(gene))
        else:
            print(f'Gene {gene} not found in the list of gene names of the GTEx data')
    # only keep the genes of interest in the dataframe
    gtex_df = gtex_df.iloc[indices_of_interest]
    # make descriptions the index
    gtex_df = gtex_df.set_index('Description')
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

    # save the dataframes
    gtex_df.to_csv('../GTEx_data/gtex_rhythmonome_genes.csv')
    gtex_high_rin_df.to_csv('../GTEx_data/gtex_rhythmonome_genes_high_RIN.csv')

    # plot the corner plots for all genes in the rhythmonome
    nSubsets = 5
    lenSubset = int(np.ceil(len(genes_of_interest) / nSubsets))
    lenSubset_last = len(genes_of_interest) - lenSubset * (nSubsets - 1)
    for iSubset in range(nSubsets):
        if iSubset == nSubsets - 1:
            subset = genes_of_interest[iSubset * lenSubset:]
        else:
            subset = genes_of_interest[iSubset * lenSubset:(iSubset + 1) * lenSubset]
            print(subset)
        subset_df = gtex_df.loc[subset]
        #  drop colum Name
        subset_df = subset_df.drop(columns=['Name'])
        subset_df = subset_df.astype(float)
        subset_df = subset_df.T
        # no extract the similar subset from the high RIN dataframe
        subset_df_high_rin = gtex_high_rin_df.loc[subset]
        # subset_df_high_rin = subset_df_high_rin.drop(columns=['Name'])
        subset_df_high_rin = subset_df_high_rin.astype(float)
        subset_df_high_rin = subset_df_high_rin.T
        # plot the corner plots for two subsets on the same figure
        fig = corner.corner(subset_df, labels=subset, color='orange', hist_kwargs={'color': 'orange'},show_contours=False,show_datapoints=True)
        corner.corner(subset_df_high_rin, labels=subset, fig=fig, color='purple', hist_kwargs={'color': 'purple'},show_contours=False,show_datapoints=True)
        # add a legend
        blue_line = mlines.Line2D([], [], color='orange', label='All RIN')
        red_line = mlines.Line2D([], [], color='purple', label='High RIN')
        # pd.plotting.scatter_matrix(subset_df, figsize=(15, 15), marker='o', color='darkorange',
        #                            hist_kwds={'bins': 20,'color': 'darkorange'}, s=15, alpha=.5)
        # #plot the high rin dat frame scatter matrix in the same plot
        # #  get current axes
        # ax = plt.gca()
        # pd.plotting.scatter_matrix(subset_df_high_rin,  marker='o', color='purple',
        #                            hist_kwds={'bins': 20,'color': 'purple'}, s=15, alpha=.5,ax=ax)
        # blue_line = mlines.Line2D([], [], color='darkorange', label='All RIN')
        # red_line = mlines.Line2D([], [], color='purple', label='High RIN')
        plt.legend(handles=[blue_line, red_line])
        plt.tight_layout()
        plt.savefig(f'Figures/GTEx_corner_plot_{iSubset}.png')
        plt.close()

    print('pause here')
