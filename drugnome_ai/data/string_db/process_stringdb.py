import pandas as pd
import numpy as np

protein_links = pd.read_csv("9606.protein.links.v11.0.txt.gz", sep=' ')
physical_links = pd.read_csv("9606.protein.physical.links.v11.0.txt.gz")
hgnc = pd.read_csv("ensp_to_hgnc.txt", sep='\t')
hgnc = hgnc[hgnc['Protein stable ID'].notnull()]
protein_names = dict(zip(hgnc['Protein stable ID'].tolist(), hgnc['Gene name'].tolist()))
datasets = {'protein_links': protein_links,
            'physical_links': physical_links}

for name, df in datasets.items():
    df[['protein1', 'protein2']] = df[['protein1', 'protein2']].apply(lambda x: x.str.lstrip('9606.'))
    df[['protein1', 'protein2']] = df[['protein1', 'protein2']].apply(lambda x: x.str.strip())

    df['gene1'] = df['protein1'].apply(lambda x: protein_names[x] if x in protein_names.keys() else None)
    df['gene2'] = df['protein2'].apply(lambda x: protein_names[x] if x in protein_names.keys() else None)

    df = df[df.gene1.notnull()]
    df = df[df.gene2.notnull()]

    # df = df[protein_links.combined_score >= 317]
    interactions = pd.DataFrame(df.gene1.unique(), columns=['Gene_Name'])

    interactions.set_index('Gene_Name', inplace=True)
    interactions['Gene_Name'] = interactions.index
    interactions.sort_index(inplace=True)

    df.set_index('gene1', inplace=True)
    df['Gene_Name'] = df.index
    df = df[df.index.isin(protein_names.values())]
    df.sort_index(inplace=True)

    interactions['interacting_genes'] = df.groupby('gene1')['gene2'].apply(list)
    interactions.to_csv(name + '_interactions.tsv.gz', sep='\t', index=False, compression='gzip')
    df.to_csv("processed_" + name + ".csv.gz", index=False, compression='gzip')
