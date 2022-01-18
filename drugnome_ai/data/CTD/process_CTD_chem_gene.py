## python script
import pandas as pd
import numpy as np

print('\n> Processing CTDbase chemical_gene interactions file...')
chem_df = pd.read_csv("CTD_chem_gene_ixns.csv.gz", engine='python', index_col=False)
d = {'# ChemicalName': "ChemicalName"}
chem_df.rename(columns = d, inplace = True)
chem_df = chem_df[chem_df.Organism=='Homo sapiens']
chem_df.InteractionActions = chem_df.InteractionActions.apply(lambda x: x.split("|"))
long_chem_df = chem_df.explode("InteractionActions")
long_chem_df.to_csv('processed_CTD_chem_gene_ixns.csv.gz', index=False, compression='gzip')