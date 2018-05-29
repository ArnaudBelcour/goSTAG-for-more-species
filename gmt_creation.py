import argparse
import numpy as np
import os
import pandas as pa
import pronto
import re

from ast import literal_eval
from tqdm import *

parser = argparse.ArgumentParser(prog = "gmt_DE_gene_creation.py")
parser.add_argument("-a", "--annot", dest = "annot_file", metavar = "FILE", help = "TSV Annotation file.", required = True)
parser.add_argument("-f", "--folder", dest = "input_folder", metavar = "FILE", help = "Folder containing De genes files.", required = True)
args = parser.parse_args()

def genes_file_creation(input_folder):
    """
    Creates a gmt file containing DE gene list associated with
    each gene in these list.

    To know what are the DE gene list file, the script look at the 
    isolating word in the file name to take only the DE gene lists file.

    If you don't know the name of the gene column, this script will show you
    all the column name, you have just to look at all these name sand isolated
    the name of the gene column.
    """
    file_paths = {}
    for file_name in os.listdir(input_folder):
        file_paths[file_name] = input_folder + '/' + file_name

    df = pa.DataFrame()
            
    for file_name in file_paths:
        df_temp = pa.read_csv(file_paths[file_name], sep='\t', header=None)
        print(df_temp.columns)
        gene_column = 0
        df_temp = df_temp[[gene_column]]
        df_temp.columns = ['Gene_Name_DE']
        row = []
        file_extension = os.path.splitext(file_name)[1]
        row.append(file_name.replace(file_extension, ""))
        row.extend(df_temp['Gene_Name_DE'].tolist())
        df = df.append([row], ignore_index=True)

    df.insert(1, 'Description', 'Genes_DE')

    df.to_csv('DE_gene.gmt', sep='\t', index=False, header=False)

def gos_files_creation(annotation_file):
    """
    This script creates a gmt file (tab-separated value file), containing each GO term
    associated with all the genes that are annotated by this GO term.

    To do this, we use a genome annotation file containing all the gene associated with
    GO terms. We also use namespaces and alternative ID of GO term.
    The namespace is to look only at one namespace (this limitation comes from goSTAG).
    The alternative ID is to avoid missing GO terms because they have multiple IDs.

    You have to set manually the namespace that you are looking at.

    """
    namespace = 'molecular_function'

    go_ontology = pronto.Ontology('http://purl.obolibrary.org/obo/go/go-basic.obo')

    # For each GO terms look to the namespaces associated with them.
    go_namespaces = {}
    for go_term in go_ontology:
        go_namespaces[go_term.id] = go_term.other['namespace'][0]

    # For each GO terms look if there is an alternative ID fo them.
    go_alt_ids = {}
    for go_term in go_ontology:
        if 'alt_id' in go_term.other:
            for go_alt in go_term.other['alt_id']:
                go_alt_ids[go_alt] = go_term.id

    # Genome file with genes associated with GO terms.
    df = pa.read_csv(annotation_file, sep='\t', header=None)
    df.columns = ['Gene_Name', 'GOs']
    df.replace(np.nan, '', inplace=True)

    gos_in_df = []
    for gos in df['GOs']:
        for go in gos.split(','):
            if go not in gos_in_df:
                gos_in_df.append(go)

    df.set_index('Gene_Name', inplace=True)

    gene_gos = []
    for gene, row in df.iterrows():
        for go in row['GOs'].split(','):
            gene_gos.append((go, gene))

    dic_go_genes = {}
    for go in tqdm(gos_in_df):
        genes = []
        for gene_go in gene_gos:
            if go != '' and go not in go_namespaces:
                go = go_alt_ids[go]
            if gene_go[0] == go and go != '' and go_namespaces[go] == namespace:
                genes.append(gene_go[1])
        if go != '':
            dic_go_genes[go] = genes

    print(len(dic_go_genes))

    delete_keys = []
    for go in dic_go_genes:
        if len(dic_go_genes[go]) < 4:
            delete_keys.append(go)

    for key in delete_keys:
        del dic_go_genes[key]
    print(len(dic_go_genes))

    df_go = pa.DataFrame.from_dict(dic_go_genes, orient='index')
    df_go.insert(0, 'Description', 'GO_terms')

    df_go.to_csv('go_gene.gmt', sep='\t', header=False)

    df.reset_index(inplace=True)
    df_query_go = pa.concat([pa.Series(row['Gene_Name'], row['GOs'].split(','))
                        for _, row in df.iterrows()]).reset_index()
    df_query_go.columns = ['GOs', 'Gene_Name']
    df_query_go = df_query_go[['Gene_Name', 'GOs']]
    df_query_go.to_csv('query_go.tsv', sep='\t', index=False)

def main(args):
    input_folder = args.input_folder
    annotation_file = args.annot_file
    genes_file_creation(input_folder)
    gos_files_creation(annotation_file)

if __name__ == '__main__':
    main(args)
