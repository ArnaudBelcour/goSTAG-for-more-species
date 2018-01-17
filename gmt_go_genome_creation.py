import numpy as np
import pandas as pa
import pronto

from ast import literal_eval
from tqdm import *

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
df = pa.read_csv('genome_file_completed.tsv', sep='\t')
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
df_go = df_go[df_go[0].str.startswith('Pldbra') == True]

df_go.to_csv('go_gene.gmt', sep='\t', header=False)
