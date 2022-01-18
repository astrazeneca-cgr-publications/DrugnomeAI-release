import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import sys
import re
import ast
from pathlib import Path

from drugnome_ai.config_class import Config


class ProcessLabels:

    def __init__(self, cfg):
        self.cfg = cfg

    def process_druggability_labels(self):
        '''
        Extract druggability tier labels
        :return: drug_labels_df
        '''
        print("\n>> Compiling druggability labels...")

        drug_labels_df = pd.read_csv(self.cfg.data_dir / 'labels/gene_druggable_labels.csv')

        drug_labels_df = drug_labels_df.dropna(subset=['Gene_Name'])
        drug_labels_df = drug_labels_df[['druggability_tier', 'Gene_Name']]
        if hasattr(self.cfg, 'tier_tag'):
            for tier in self.cfg.tier_tag:
                print(tier)
                drug_labels_df.replace({str(tier): 1}, inplace=True)

            if len(self.cfg.tier_tag) != 4:
                for word in ['Tier 1', 'Tier 2', 'Tier 3A', 'Tier 3B']:
                    if word not in self.cfg.tier_tag:
                        drug_labels_df.replace({str(word): 0}, inplace=True)

        drug_labels_df['known_gene'] = drug_labels_df['druggability_tier']
        drug_labels_df.drop(columns = 'druggability_tier', inplace = True)
        return drug_labels_df


    def process_pharos_labels(self):
        '''
        Extract druggability related labels derived from Pharos
        :return: pharos_df
        '''
        print("\n>> Compiling Pharos labels...")

        if not hasattr(self.cfg, 'config_file'):
            pharos_df = pd.read_csv(self.cfg.data_dir / 'PHAROS/pharos_GF_wINDEX.csv')

            pharos_df = pharos_df[['Gene_Name', 'idgTDL']]
            pharos_df['Gene_Name'].drop_duplicates(inplace=True)

            for tier in self.cfg.pharos_tag:
                print(tier)
                pharos_df.replace({str(tier): 1}, inplace=True)

            if len(self.cfg.pharos_tag) != 4:
                for word in ['Tclin', 'Tchem', 'Tbio', 'Tdark']:
                    if word not in self.cfg.pharos_tag:
                        pharos_df.replace({str(word): 0}, inplace=True)

        else:
            pharos_df = pd.read_csv(self.cfg.processed_label_table, sep='\t', index_col=None)

        if 'idgTDL' in pharos_df.columns:
            pharos_df['known_gene'] = pharos_df['idgTDL']
            pharos_df.drop(columns='idgTDL', inplace=True)

        return pharos_df


    def process_hgnc_names(self):
        '''
        Extract gene names derived from HUGO Gene Nomenclature Committee resource
        :return: gene_df
        '''
        gene_df = pd.read_csv(self.cfg.data_dir / 'PHAROS/pharos_GF_wINDEX.csv', index_col=None)

        gene_df = gene_df[['Gene_Name']]
        gene_df.dropna(inplace=True)
        gene_df.drop_duplicates(inplace=True)

        return gene_df

############ NEW
    def provide_seed_genes(self):

        pos_labs_df = pd.read_csv(self.cfg.processed_label_table, sep='\t', index_col=None)
        print(pos_labs_df.head())

        pos_labs_df = pos_labs_df.loc[pos_labs_df['known_gene'] == 1]

        return pos_labs_df['Gene_Name'].tolist()


    def run_all(self):

        gene_names = self.process_hgnc_names()

        if hasattr(self.cfg, 'pharos_tag'):
            drug_labels_df = self.process_pharos_labels()
            print(drug_labels_df.columns)

        elif hasattr(self.cfg, 'tier_tag'):
            drug_labels_df = self.process_druggability_labels()
            print(drug_labels_df.columns)

        else:
            try:
                drug_labels_df = pd.read_csv(self.cfg.custom_known_genes_file, header=None)
                drug_labels_df.columns = ['Gene_Name']
                drug_labels_df['known_gene'] = 1
                print(drug_labels_df.columns)
            except:
                sys.exit(
                    "[Error] Could not read input file with custom known genes list.\nPlease provide a file with HGNC gene names (each in a separate line) with the -k option when calling drugnomeai or provide label type instead with -t or -p option.")

        drug_labels_df['known_gene'].fillna(0, inplace=True)
        labels_df = pd.merge(gene_names, drug_labels_df, how='left', on='Gene_Name')
        labels_df.drop_duplicates(inplace=True)
        labels_df['known_gene'].fillna(0, inplace=True)

        if self.cfg.random_seeds:
            print(labels_df.loc[labels_df['known_gene'] == 1, :].shape)

            total_seed_genes = labels_df.loc[ labels_df['known_gene'] == 1, :].shape[0]
            # reset known genes labels with '0' value for all genes
            labels_df.loc[:, 'known_gene'] = 0
            print(labels_df.loc[ labels_df['known_gene'] == 1, :].shape)

            # select random indexes
            random_seed_indexes = np.random.choice(list(range(labels_df.shape[0])), size=total_seed_genes, replace=False).tolist()
            print(len(random_seed_indexes))

            # assign '1' value for 'known_gene' label to random genes indicated by the generated random indexes
            labels_df.loc[ random_seed_indexes, 'known_gene'] = 1
            print(labels_df.loc[ labels_df['known_gene'] == 1, :].shape)


        labels_df.to_csv(self.cfg.processed_label_table, sep='\t', index=None)
        print("Saved to {0}".format(self.cfg.processed_label_table))

        print(labels_df.shape)

if __name__ == '__main__':

    out_dir = sys.argv[1]
    cfg = Config(out_dir)
    cfg.tier_tag = ['Tier 1']
    proc = ProcessLabels(cfg)
    proc.run_all()

