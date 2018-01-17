import os
import pandas as pa
import re

"""
This file creates a gmt file containing DE gene list associated with
each gene in these list.

To know what are the DE gene list file, the script look at the 
isolating word in the file name to take only the DE gene lists file.

If you don't know the name of the gene column, this script will show you
all the column name, you have just to look at all these name sand isolated
the name of the gene column.
"""

isolating_word = 'diff'

file_names = []
for file_name in os.listdir():
    if isolating_word in file_name:
        file_names.append(file_name)

df = pa.DataFrame()
        
for file_name in file_names:
    df_temp = pa.read_csv(file_name, sep='\t')
    print(df_temp.columns)
    gene_column = 'Unnamed: 0'
    df_temp = df_temp[[gene_column]]
    df_temp.columns = ['Gene_Name_DE']
    row = []
    file_extension = os.path.splitext(file_name)[1]
    row.append(file_name.replace(file_extension, ""))
    row.extend(df_temp['Gene_Name_DE'].tolist())
    df = df.append([row], ignore_index=True)

df.insert(1, 'Description', 'Genes_DE')

df[0] = df[0].str.replace(isolating_word, '')
df.to_csv('DE_gene.gmt', sep='\t', index=False, header=False)
