
from setup import *
import sympy as sym
from scipy.io import mmread
import matrixconverters as mcs

def get_gene_indeces(gene_names, genes):
    """
    Get the indeces of the gene in the genes dataframe
    :param gene_names: list of str, names of the gene (capitalised)
    :param genes: pd.DataFrame, dataframe containing all gene names
    :return: np.array, indeces of the genes of interest
    """
    # if gene_names are not capitalised, capitalise them
    gene_names = [gene.upper() for gene in gene_names]
    gene_indeces = []
    for gene in gene_names:
        gene_indeces.append(genes[genes['gene_name'].str.contains(gene)].index[0])
    gene_indeces = np.array(gene_indeces)
    return gene_indeces


if __name__ == '__main__':
    ####################################################################################################################
    # we only needed to run this to undestand the data format and how it matches other tables. We can skip this now
    # filename = '../SingleCellGene_data/gene_sorted-matrix.mtx'
    # try_reading_data = mmread(filename)
    # # get number of columns and rows
    # n_rows = try_reading_data.shape[0]
    # n_columns = try_reading_data.shape[1]
    # print('Number of rows:', n_rows)
    # print('Number of columns:', n_columns)
    # del try_reading_data
    ####################################################################################################################
    # get gene names
    filename = '../SingleCellGene_data/genes_v2.tsv'
    with open(filename, 'r') as f:
        genes = f.readlines()
    genes = [gene.strip().split('\t') for gene in genes]
    genes = pd.DataFrame(genes, columns=['gene_id', 'gene_name'])
    genes.to_csv('../SingleCellGene_data/genes_names_readable.csv', index=False)
    print(f'Number of gene ids:{genes.shape[0]}')
    unique_gene_names = genes['gene_name'].unique()
    n_genes = len(unique_gene_names)
    print(f'Number of unique gene names:{n_genes}')
    # check how many gene names contain 'KCNH2'
    gene_name = 'KCNH2'
    gene_names_with_herg = genes[genes['gene_name'].str.contains(gene_name)]
    # read the metadata matrix
    filename = '../SingleCellGene_data/meta.data.v3.txt'
    with open(filename, 'r') as f:
        metadata = f.readlines()
    # convert the metadata to pandas dataframe
    metadata = [meta.strip().split('\t') for meta in metadata]
    metadata = pd.DataFrame(metadata[1:], columns=metadata[0])
    # drop the first row of metadata
    metadata = metadata.drop(0, axis=0).reset_index(drop=True)
    # get the length of metadata
    n_metadata = metadata.shape[0]
    # save the metadata to csv file
    metadata.to_csv('../SingleCellGene_data/metadata_readable.csv', index=False)
    # get the number of unique cell names from metadata
    cell_names = metadata['NAME'].unique()
    n_cells = len(cell_names)
    print(f'Number of cells:{n_cells}')
    ####################################################################################################################
    # we now know that the number of rows in the matrix is equal to the number of genes
    # the number of columns in the matrix is equal to the number of cells

    # from metadata we are only interested in cardiomyocytes, so we can get those indeces
    cardiomyocyte_cell_index = metadata[metadata['Cluster'].str.contains('Cardiomyocyte')].index
    cluster_names_to_store = metadata['Cluster'].iloc[cardiomyocyte_cell_index].tolist()
    cell_names_to_store = metadata['NAME'].iloc[cardiomyocyte_cell_index].tolist()
    # maximum value of cardiomyocyte_cell_index
    print(f'Maximum value of cardiomyocyte_cell_index:{max(cardiomyocyte_cell_index)}')
    # ^^ note that these contain cytoplasmic and atrial cardiomyocytes as well
    # get all the genes that were in the rhythmome (plus RYR2 and LAMP1 as negative controls)
    genes_of_interest = genes_in_rhythmonome + ['LAMP1']
    gene_indeces = get_gene_indeces(genes_of_interest, genes)
    # order genes of interest based on the order of gene_indeces (to make matrix extraction easier)
    gene_indeces = np.sort(gene_indeces)
    genes_of_interest = genes['gene_name'].iloc[gene_indeces].tolist()
    # maximum value of gene_indeces
    print(f'Maximum value of gene_indeces:{max(gene_indeces)}')
    # check for missing gene names in the list of genes
    missing_genes = [gene for gene in genes_in_rhythmonome if gene not in genes_of_interest]
    if missing_genes is not None:
        print(f'Missing genes in single cell data: {missing_genes}')
    ####################################################################################################################
    # now read the big matrix, and only extract columns corresponding to cardiomyocytes and rows corresponding to genes of interest
    # read the matrix
    filename = '../SingleCellGene_data/gene_sorted-matrix.mtx'
    data = mmread(filename)
    # get the number of columns and the number of rows
    n_rows = data.shape[0]
    n_columns = data.shape[1]
    print(f'Number of rows in the matrix:{n_rows}')
    print(f'Number of columns in the matrix:{n_columns}')
    # check which shape exceeds the number cardiomyocyte_cell_index maximum value
    # extract the data for genes only, then for cardiomyocytes only (weird indexing issue if we try do do both)
    data = data.tocsr()[gene_indeces, :]
    data = data[:, cardiomyocyte_cell_index]
    # convert the daa to pandas dataframe, make sure that indeces are cells and columns are genes
    data_df = pd.DataFrame(data.toarray().transpose())
    # add the gene names to the dataframe
    data_df.columns = genes_of_interest
    # add the cell names to the dataframe
    data_df['Cell_name'] = cell_names_to_store
    data_df['Cluster'] = cluster_names_to_store
    # move the cell name in cluster to the first two columns
    cols = data_df.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    data_df = data_df[cols]
    # save the dataframe
    data_df.to_csv('../SingleCellGene_data/cardiomyocytes_rhythmonome_only.csv', index=False)

    print('pause here')