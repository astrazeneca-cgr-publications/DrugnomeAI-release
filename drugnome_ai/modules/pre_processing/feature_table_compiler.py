import sys
import re
import numpy as np
import pandas as pd
import math
from functools import reduce
from drugnome_ai.modules.pre_processing.data_compilation import (process_generic_features,
                                                               process_druggability_features,
                                                               process_labels,
                                                               process_features_filtered_by_disease
                                                               )

from drugnome_ai.config_class import Config



class FeatureTableCompiler:

    def __init__(self, cfg):
        self.cfg = cfg
        self.full_features_df = None
        print("Feature Table Compiler launched")


    def compile_feature_tables_per_class(self):
        print("\n>> Compiling feature tables per class...")

        # filtered by disease
        proc = process_features_filtered_by_disease.ProcessFeaturesFilteredByDisease(self.cfg)
        proc.run_all()

        # labels
        labeler = process_labels.ProcessLabels(self.cfg)
        labeler.run_all()
        
        # generic
        proc = process_generic_features.ProcessGenericFeatures(self.cfg)
        proc.run_all()

        # druggability
        proc = process_druggability_features.ProcessDruggabilityFeatures(self.cfg)
        proc.run_all(labeler)

    def combine_all_feature_tables(self):
        print("\n>> Combining all feature tables together...")

        data_sources = {
            'pharos': self.cfg.druggability_feature_table,
            'inter_pro' : self.cfg.druggability_feature_table
        }

        datasets = []
        datasets.append(pd.read_csv(self.cfg.processed_label_table, sep='\t', index_col=None))

        if hasattr(self.cfg, 'genic_intol_option') or hasattr(self.cfg, 'mantis_data_option'):
            datasets.append(pd.read_csv(self.cfg.generic_feature_table, sep='\t', index_col=None))

        for source in self.cfg.data_source:
            path = data_sources[source]
            if ('pharos' in self.cfg.data_source) & (source == 'inter_pro'):
                print('InterPro already loaded')
                pass
            else:
                datasets.append(pd.read_csv(path, sep='\t'))

        generic_features_df = reduce(lambda df1, df2: pd.merge(df1, df2, how='outer', on='Gene_Name'), datasets)

        print(generic_features_df.shape)

        filtered_by_tissue_df = pd.read_csv(self.cfg.filtered_by_disease_feature_table, sep='\t')
        self.full_features_df = pd.merge(generic_features_df, filtered_by_tissue_df, how='left', on='Gene_Name')

	    # dropping missing gene names
        self.full_features_df.drop_duplicates(subset='Gene_Name', inplace=True)
        self.full_features_df.reset_index(drop=True, inplace=True)

        if not 'known_gene' in self.full_features_df:
            print("Not known gene column found")

        if 'pharos' in self.cfg.data_source:
            self.full_features_df.drop(columns = 'pharos_index', inplace = True)
            self.full_features_df[['Re_L1_seed_genes_overlap', 'Re_L2_seed_genes_overlap']].fillna(0, inplace=True)

        self.full_features_df['known_gene'].fillna(0, inplace=True)

    def inspect_missing_data(self, df, verbose=False):

        # missing data
        total = df.isnull().sum().sort_values(ascending=False)
        percent = (df.isnull().sum() / df.isnull().count()).sort_values(ascending=False)

        missing_data = pd.concat([total, percent], axis=1, keys=['Total', 'Percent'])
        missing_data = missing_data.loc[missing_data['Percent'] > 0]
        missing_data['Percent'] = missing_data['Percent'].apply(lambda x: round(x * 100, 2))
        print("Number of features with missing data: {0}".format(missing_data.shape[0]))

        full_missing_data_df = missing_data.copy()

        if self.cfg.create_plots and missing_data.shape[0] > 0:

            if self.cfg.generic_classifier:
                # collapse GO features into a single feature
                indexes_to_collapse = [c for c in missing_data.index.values if c.startswith('GO_')]
                collapsed_row = missing_data.loc[indexes_to_collapse[0]].copy()
                collapsed_row.name = 'GO-collapsed_features'
                missing_data = missing_data.append(collapsed_row)

                missing_data.drop(indexes_to_collapse, axis=0, inplace=True)

                # collapse ProteinAtlas features into a single feature
                indexes_to_collapse = [c for c in missing_data.index.values if c.startswith('ProteinAtlas_')]
                collapsed_row = missing_data.loc[indexes_to_collapse[0]].copy()
                collapsed_row.name = 'ProteinAtlas-collapsed_features'
                missing_data = missing_data.append(collapsed_row)

                missing_data.drop(indexes_to_collapse, axis=0, inplace=True)

                # collapse GTEx features into a single feature
                indexes_to_collapse = [c for c in missing_data.index.values if c.startswith('GTEx_')]
                collapsed_row = missing_data.loc[indexes_to_collapse[0]].copy()
                collapsed_row.name = 'GTEx-collapsed_features'
                missing_data = missing_data.append(collapsed_row)

                missing_data.drop(indexes_to_collapse, axis=0, inplace=True)

                missing_data.sort_values(by='Percent', ascending=False, inplace=True)


            ax = missing_data.reset_index().plot.barh(x='index', y='Percent',
                                                      align='center', color='#4292c6', ecolor='black',
                                                      fontsize=8, figsize=(10, 15))

            vline_thres = self.cfg.missing_data_thres * 100
            ax.axvline(vline_thres, color="#de2d26", linestyle='-.', linewidth=0.6)
            ax.invert_yaxis()  # labels read top-to-bottom
            ax.set_xlabel('Missing data % ratio')
            ax.xaxis.set_label_position('top')
            ax.xaxis.tick_top()

            fig = ax.get_figure()
            plot_filepath = str(self.cfg.eda_out / 'missing_data_ratios.pdf')
            fig.savefig(plot_filepath, format='pdf', bbox_inches='tight')

        if verbose:
            missing_data_str = "|\tFeature\t|\tTotal\t|\tPercent\t|\n"
            for index, row in missing_data.iterrows():
                missing_data_str += index + "\t" + str(row['Total']) + "\t" + str(row['Percent']) + "\t\n"
            print(missing_data_str)

        missing_data = full_missing_data_df

        return missing_data


    def drop_features_w_missing_data(self, df, missing_data, missing_data_thres=0.99):
        '''
        Drop features with high ratio of missing data

        :param df: 
        :param missing_data: 
        :param missing_data_thres: 
        :return: 
        '''
        print('\n>> Removing features with high ratio of missing data...')
        print('\n>> Removing unkown genes')
        

        self.full_features_df['Gene_Name'] = self.full_features_df['Gene_Name'].dropna()

        missing_data_thres *= 100
        missing_data_elements = missing_data.loc[missing_data['Percent'] > missing_data_thres].index.values

        elems_to_drop = ', '.join(missing_data_elements)

        df = df.drop(missing_data_elements, axis=1)
        print('Dropped {0} features with more than {1}% missing values'.format(str(len(missing_data_elements)),
                                                                               str(missing_data_thres)))
        print(elems_to_drop)

        return (df)


    def impute_nas_w_zeros(self, df):

        valid_feature_substrings = ['GO_', 'ProteinAtlas_', 'DGIdb_', 'GTEx_',
                                    'GWAS_', 'IPR_', 'Re_', 'in_web', 'string_db', 'ctd', 'OMIM']
        go_cols = [col for col in df.columns if any(s in col for s in valid_feature_substrings)]

        replace_with_zero_elements = ['known_gene', 'glom_FDR', 'glom_Pr_of_no_eQTL', 'glom_Exp_num_of_eQTLs',
                                      'tub_FDR', 'tub_Pr_of_no_eQTL', 'tub_Exp_num_of_eQTLs', 'MGI_mouse_knockout_feature',
                                      'GOA_Kidney_Research_Priority', 'DAPPLE_perc_core_overlap', 'Inferred_perc_core_overlap',
                                      'Experimental_perc_core_overlap', 'MGI_essential_gene', 'ProteinAtlas_gene_expr_levels',
                                      'tubular_expr_flag', 'glomerular_expr_flag', 'CKDdb_num_of_studies', 'CKDdb_Disease',
                                      'ProteinAtlas_RNA_expression_TMP', 'essential_mouse_knockout', 'non_essential_mouse_knockout',
                                      'platelets_eQTL', 'HT_eQTL_hits', 'CAD_eQTL_hits', 'adipose_cis_eQTL', 'adipose_GWAS_locus',
                                      'druggability_tier', 'drug_claim_name', 'MouseGenes_Non-essential',
                                      'MouseGenes_Essential', 'MouseGenes_Overexpressed_in_the_brain', 'ppiCount'
                                     ]

        replace_with_zero_elements.extend(go_cols)

        for col in replace_with_zero_elements:
            if col in df.columns:
                df[col].fillna(0, inplace=True)

        if 'ProteinAtlas_gene_expr_levels' in df.columns:
            df.loc[ df.ProteinAtlas_gene_expr_levels.isin([0, '0']), 'ProteinAtlas_gene_expr_levels'] = 'Not_detected'

        return df


    def impute_nas_w_median(self, df):
        valid_feature_substrings = ['GnomAD_', 'GenicIntolerance_']
        gnomad_cols = [col for col in df.columns if any(s in col for s in valid_feature_substrings)]

        replace_with_median_elements = ['ExAC_dup.score', 'ExAC_del.score', 'ExAC_dup.sing.score', 'ExAC_del.sing.score',
                                        'ExAC_dup.sing', 'ExAC_del.sing', 'ExAC_num_targ', 'ExAC_dup', 'ExAC_del',
                                        'ExAC_mean_rd', 'ExAC_gc_content', 'ExAC_complexity', 'ExAC_cds_len',
                                        'ExAC_gene_length', 'ExAC_flag', 'ExAC_segdups', 'ExAC_dip', 'ExAC_cnv.score',
                                        'mut_prob_splice_site', 'GeneSize',
                                        'hpa_prot_spec', 'gene_length', 'gc_content', 'complexity', 'cds_len', 'num_targ',
                                        'cnv.score', 'mean_rd', 'hpm_gene_spec', 'hpm_prot_spec', 'MouseGenes_GC%', 'novelty', 'hpa_RNA_spec',
                                        'antibodyCount', 'monoclonalCount', 'uniprot_seq_len']

        replace_with_median_elements.extend(gnomad_cols)

        for col in replace_with_median_elements:
            if col in df.columns:
                df[col].fillna(0, inplace=True)

        return df


    def impute_missing_data(self):

        print("\n>> Imputing missing data...")
        self.impute_nas_w_zeros(self.full_features_df)
        self.impute_nas_w_median(self.full_features_df)
        self.full_features_df.dropna(subset = ['Gene_Name'], inplace=True)



    def verify_no_missing_data(self, df):

        missing_data = df.isnull().sum().sort_values(ascending=False)
        missing_data = missing_data[ missing_data != 0]
        print("Number of features with missing data: {0}".format(missing_data.shape[0]))
        if missing_data.shape[0] != 0:
            print(missing_data)
            print("[Error]: Feature table contains missing data. Aborting...")
            sys.exit()
        else:
            print("All missing data have been successfully imputed.")


    def run(self):
        # compile feature tables per class (generic/filtered-by-disease/disease-specific)
        self.compile_feature_tables_per_class()

        # merge all tables
        self.combine_all_feature_tables()

        # check for missing data
        missing_data = self.inspect_missing_data(self.full_features_df, verbose=True)
        if self.cfg.drop_missing_data_features:
            self.full_features_df = self.drop_features_w_missing_data(self.full_features_df, missing_data, missing_data_thres=self.cfg.missing_data_thres)

        # impute missing data
        self.impute_missing_data()
        self.verify_no_missing_data(self.full_features_df)

        # move 'known_gene' column to end
        tmp_known_gene = self.full_features_df[self.cfg.Y]
        self.full_features_df.drop(self.cfg.Y, axis=1, inplace=True)
        self.full_features_df[self.cfg.Y] = tmp_known_gene
        print(self.full_features_df.shape)

        self.full_features_df.to_csv(self.cfg.complete_feature_table, sep='\t', index=None)
        print("Saved full feature table (after imputation) to {0}".format(self.cfg.complete_feature_table))


if __name__ == '__main__':

    config_file = sys.argv[1] #'../../../config.yaml'
    out_dir = sys.argv[2]
    cfg = Config(out_dir, config_file)
    cfg.tier_tag = ['Tier 1']
    feat_compiler = FeatureTableCompiler(cfg)
    feat_compiler.run()

