import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import sys
import re
import ast
from functools import reduce
from scipy.stats import hmean
from pathlib import Path
from sklearn.preprocessing import MultiLabelBinarizer

from drugnome_ai.config_class import Config


class ProcessDruggabilityFeatures:

    def __init__(self, cfg):
        self.cfg = cfg

    def process_pharos(self, save_to_file=False):
        '''
        Process PHAROS features pulled from API (June 2020 - version 5.4.0)
        :return: pharos_df
        '''
        print("\n>> Compiling PHAROS...")

        pharos_df = pd.read_csv(self.cfg.data_dir / 'PHAROS/pharos_GF_wINDEX.csv')

        pharos_df = pharos_df.drop(columns='AccessionKeys')

        if not hasattr(self.cfg, 'extra_features_option'):
            pharos_df = pharos_df.drop(columns='novelty')

        pharos_df[['ppiCount']] = pharos_df[['ppiCount']].fillna(0)

        pharos_df['Gene_Name'].dropna(inplace=True)
        pharos_df['Gene_Name'].drop_duplicates(inplace=True)

        # if hasattr(self.cfg, 'pharos_tag'):
        #     for tier in self.cfg.pharos_tag:
        #         print(tier)
        #         pharos_df.replace({str(tier): 1}, inplace=True)
        #
        #     if len(self.cfg.pharos_tag) != 4:
        #         for word in ['Tclin', 'Tchem', 'Tbio', 'Tdark']:
        #             if word not in self.cfg.pharos_tag:
        #                 pharos_df.replace({str(word): 0}, inplace=True)
        # else:
        pharos_df.drop(['idgTDL'], axis=1, inplace=True)

        if save_to_file:
            pharos_df.to_csv(self.cfg.data_dir / 'PHAROS/compiled_pharos_features.tsv', sep='\t', index=None)

        return pharos_df

    def process_dgidb(self, save_to_file=False):
        '''
        Process features extracted from drug-gene interaction database (https://www.dgidb.org/)
        :return: dgidb_count_df
        '''
        print("\n>> Compiling DGIdb features...")

        dgidb_df = pd.read_csv(self.cfg.data_dir / 'DGIdb/interactions.csv')

        dgidb_df.rename(columns= {'gene_name': 'Gene_Name'}, inplace=True)
        dgidb_count_df = dgidb_df.groupby("Gene_Name").count()
        dgidb_count_df = dgidb_count_df['interaction_types']

        dgidb_count_df.rename('DGIdb_interaction_types', inplace=True)

        print(dgidb_count_df.shape)
        if save_to_file:
            dgidb_count_df.to_csv(self.cfg.data_dir / 'DGIdb/dgidb_compiled_feature_table.tsv', sep='\t', index=False)

        return dgidb_count_df

    def process_reactome(self, seed_genes):
        '''
        Process features related to protein-protein interactions (PPIs)
        between genes being level-1 and level-2 neighbours
        :return: reactome_ppi_df
        '''

        print("\n>> Compiling Reactome features...")
        seed_genes = set(seed_genes)

        reactome_ppi_df = pd.read_csv(self.cfg.data_dir / 'Reactome/reactome_ppi_final.txt', sep='\t',
                                      index_col=None)
        reactome_ppi_df.rename(columns={'Interactions': 'interacting_genes'}, inplace=True)

        #######
        tmp_df = reactome_ppi_df[:]
        reactome_genes = reactome_ppi_df.interacting_genes.apply(lambda x: list(ast.literal_eval(x)))
        reactome_genes = set(sum(reactome_genes, []))
        index_genes = reactome_ppi_df.Gene_Name.unique().tolist()
        new_genes = list(set(reactome_genes) - set(index_genes))
        new_genes_dict = dict(zip(new_genes, [None] * len(new_genes)))

        gene_present = (lambda row, gene: gene in eval(row))

        for gene in new_genes:
            tmp_df['check'] = tmp_df.interacting_genes.apply(lambda x: gene_present(x, gene))
            new_genes_dict[gene] = str(tmp_df.Gene_Name[tmp_df.check == True].tolist())

        new_genes_df = pd.DataFrame(new_genes_dict.items(), columns=['Gene_Name', 'interacting_genes'])
        reactome_ppi_df = pd.concat([reactome_ppi_df, new_genes_df])
        ######

        reactome_ppi_df['Re_L1_seed_genes_overlap'] = reactome_ppi_df['interacting_genes'].apply(
            lambda row: self.get_seed_genes_overlap(row, seed_genes))
        print('>> 1st neighbour overlap calculations complete.')

        reactome_ppi_df['Re_L2_seed_genes_overlap'] = self.layer_two_overlap(reactome_ppi_df, seed_genes)
        print('>> 2nd neighbour overlap calculations complete.')

        reactome_ppi_df.drop(['interacting_genes'], axis=1, inplace=True)
        print('\n\nreactome_ppi_df:', reactome_ppi_df.head())
        print(reactome_ppi_df.shape)

        reactome_ppi_df.fillna(0, inplace=True)
        return reactome_ppi_df

    # Generic functions for calculating features from PPI networks
    def get_seed_genes_overlap(self, interacting_genes, seed_genes):
        interacting_genes = eval(interacting_genes)
        overlapping_genes = list(set(list(interacting_genes)) & seed_genes)
        perc_overlap = len(overlapping_genes) / len(interacting_genes)
        return perc_overlap

    def layer_two_overlap(self, interacting_df, seed_genes):

        interacting_df.index = interacting_df['Gene_Name']
        interacting_df.index.name = 'gene_index'
        print(interacting_df.head())

        interacting_df['interacting_genes'] = interacting_df['interacting_genes'].apply(
            lambda x: set(ast.literal_eval(x)))

        total_store = []

        for interactions in interacting_df['interacting_genes']:

            tot_interacting_genes = 0
            tot_overlapping_genes = 0

            for gene in interactions:
                layer_two_interactions = interacting_df.loc[gene, 'interacting_genes']
                # print(layer_two_interactions)

                tot_interacting_genes += len(layer_two_interactions)
                # print('\ntot_interacting_genes:', tot_interacting_genes)

                tot_overlapping_genes += len(layer_two_interactions & seed_genes)
            # print('tot_overlapping_genes:', tot_overlapping_genes)

            layer_two_overlap = tot_overlapping_genes / tot_interacting_genes
            total_store.append(layer_two_overlap)

        return total_store

    def compile_interaction_table(self, file, source, prefix, seed_genes):
        tmp_df = pd.read_csv(str(self.cfg.data_dir) + '/' + source + '/' + file, sep='\t', index_col=None)
        print(tmp_df.head())
        tmp_df[source + '_L1_' + prefix + '_seed_genes_overlap'] = tmp_df['interacting_genes'].apply(
            lambda row: self.get_seed_genes_overlap(row, seed_genes))
        print('>> 1st neighbour overlap calucations complete.')

        tmp_df[source + '_L2_' + prefix + '_seed_genes_overlap'] = self.layer_two_overlap(tmp_df, seed_genes)
        print('>> 2nd neighbour overlap calucations complete.')

        tmp_df.drop(['interacting_genes'], axis=1, inplace=True)
        print('\n\ntmp_df:', tmp_df.head())
        print(tmp_df.shape)

        return tmp_df

    def process_inweb(self, seed_genes):
        '''
        Process features related to experimentally derived PPIs
        between genes being level-1 and level-2 neighbours
        :return: inweb_df
        '''

        print("\n>> Compiling InWeb features...")
        seed_genes = set(seed_genes)

        inferred_df = self.compile_interaction_table('inferred_pairwise_interactions.tsv', 'in_web', 'inferred',
                                                     seed_genes)
        experim_df = self.compile_interaction_table('experimental_pairwise_interactions.tsv', 'in_web', 'experimental',
                                                    seed_genes)

        inweb_df = pd.merge(experim_df, inferred_df, how='outer', left_on='Gene_Name', right_on='Gene_Name')
        inweb_df.fillna(0, inplace=True)

        print(inweb_df.head())
        print(inweb_df.tail())

        return inweb_df

    def process_interpro_domain(self):
        '''
        Process domain related features derived from combining InterPro with Pharod Generic Features
        :return: domain_df
        '''
        print("\n>> Compiling InterPro domain features...")

        # domain_df = pd.read_pickle(self.cfg.data_dir / 'interpro/domain.pkl')
        domain_df = pd.read_csv(self.cfg.data_dir / 'interpro/domain.csv')
        gene_df = pd.read_csv(self.cfg.data_dir / 'PHAROS/pharos_GF_wINDEX.csv', index_col=None)

        domain_df = domain_df[['Gene_Name', 'name']]
        gene_df = gene_df[['Gene_Name']]

        gene_df.drop_duplicates(inplace=True)

        mlb = MultiLabelBinarizer()

        # create boolean mask matched non NaNs values
        mask = domain_df['name'].notnull()

        # filter by boolean indexing
        arr = mlb.fit_transform(domain_df.pop('name').loc[mask].dropna().apply(eval))

        # create DataFrame and add missing (NaN)s index values
        multilabel_df = (pd.DataFrame(arr, index=domain_df.index[mask], columns=mlb.classes_) \
                         .reindex(domain_df.index, fill_value=0))

        domain_df = domain_df[['Gene_Name']].join(multilabel_df)

        # domain_df = domain_df.join(pd.DataFrame(mlb.fit_transform(domain_df.pop('name').dropna().str.strip('[]').str.split(',')),
        #                                         columns=mlb.classes_, index=domain_df.index))

        name_store = domain_df['Gene_Name']

        domain_df.drop('Gene_Name', inplace=True, axis=1)
        domain_df.drop([col for col, val in domain_df.sum().iteritems() if val < 97], axis=1, inplace=True)
        domain_df['Gene_Name'] = name_store

        domain_df = pd.merge(gene_df, domain_df, how='left', left_on='Gene_Name', right_on='Gene_Name')

        domain_df.drop_duplicates(inplace=True)
        domain_df.columns = ['IPR_d_' + x if x != 'Gene_Name' else x for x in domain_df.columns]
        return domain_df

    def process_interpro_family(self):
        '''
        Process family related features derived from combining InterPro with Pharod Generic Features
        :return: domain_df
        '''
        print("\n>> Compiling InterPro family features...")

        # family_df = pd.read_pickle(self.cfg.data_dir / 'interpro/family.pkl')
        family_df = pd.read_csv(self.cfg.data_dir / 'interpro/family.csv')
        gene_df = pd.read_csv(self.cfg.data_dir / 'PHAROS/pharos_GF_wINDEX.csv')

        family_df = family_df[['Gene_Name', 'name']]
        gene_df = gene_df[['Gene_Name']]

        gene_df.drop_duplicates(inplace=True)

        mlb = MultiLabelBinarizer()
        # create boolean mask matched non NaNs values
        mask = family_df['name'].notnull()

        # filter by boolean indexing
        arr = mlb.fit_transform(family_df.pop('name').loc[mask].dropna().apply(eval))

        # create DataFrame and add missing (NaN)s index values
        multilabel_df = (pd.DataFrame(arr, index=family_df.index[mask], columns=mlb.classes_) \
                         .reindex(family_df.index, fill_value=0))

        family_df = family_df[['Gene_Name']].join(multilabel_df)

        # family_df = family_df.join(pd.DataFrame(mlb.fit_transform(family_df.pop('name')),
        #                                         columns=mlb.classes_, index=family_df.index))

        name_store = family_df['Gene_Name']

        family_df.drop('Gene_Name', inplace=True, axis=1)
        family_df.drop([col for col, val in family_df.sum().iteritems() if val < 31], axis=1, inplace=True)
        family_df['Gene_Name'] = name_store

        # Val > 30 gives us the 20 most dense domains

        family_df = pd.merge(gene_df, family_df, how='left', left_on='Gene_Name', right_on='Gene_Name')

        family_df.drop_duplicates(inplace=True)
        family_df.columns = ['IPR_f_' + x if x != 'Gene_Name' else x for x in family_df.columns]

        return family_df

    def process_interpro_super_family(self):
        '''
        Process super family related features derived from InterPro with Pharod Generic Features
        :return: domain_df
        '''
        print("\n>> Compiling InterPro super family features...")

        # super_family_df = pd.read_pickle(self.cfg.data_dir / 'interpro/homologous_superfamily.pkl')
        super_family_df = pd.read_csv(self.cfg.data_dir / 'interpro/homologous_superfamily.csv')
        gene_df = pd.read_csv(self.cfg.data_dir / 'PHAROS/pharos_GF_wINDEX.csv', index_col=None)

        super_family_df = super_family_df[['Gene_Name', 'name']]
        gene_df = gene_df[['Gene_Name']]
        gene_df.drop_duplicates(inplace=True)

        mlb = MultiLabelBinarizer()

        # create boolean mask matched non NaNs values
        mask = super_family_df['name'].notnull()

        # filter by boolean indexing
        arr = mlb.fit_transform(super_family_df.pop('name').loc[mask].dropna().apply(eval))

        # create DataFrame and add missing (NaN)s index values
        multilabel_df = (pd.DataFrame(arr, index=super_family_df.index[mask], columns=mlb.classes_) \
                         .reindex(super_family_df.index, fill_value=0))

        super_family_df = super_family_df[['Gene_Name']].join(multilabel_df)

        # super_family_df = super_family_df.join(
        #     pd.DataFrame(mlb.fit_transform(super_family_df.pop('name')),
        #                  columns=mlb.classes_, index=super_family_df.index))

        name_store = super_family_df['Gene_Name']

        super_family_df.drop('Gene_Name', inplace=True, axis=1)
        super_family_df.drop([col for col, val in super_family_df.sum().iteritems() if val < 100], axis=1, inplace=True)
        super_family_df['Gene_Name'] = name_store

        # Val > 30 gives us the 20 most dense domains

        super_family_df = pd.merge(gene_df, super_family_df, how='left', left_on='Gene_Name', right_on='Gene_Name')

        ###
        super_family_df.drop_duplicates(inplace=True)
        super_family_df.columns = ['IPR_sf_' + x if x != 'Gene_Name' else x for x in super_family_df.columns]

        return super_family_df

    def process_ctd(self):

        def select_features(data, grouping_col, threshold):  # 300 for pathways, 500 for chemicals
            distribution = data.groupby([grouping_col]).GeneID.nunique().to_frame()
            #         distribution['share'] = distribution.GeneID/data.GeneID.nunique()
            distribution.sort_values('GeneID', ascending=False, inplace=True)
            selected_values_as_features = distribution[distribution.GeneID >= distribution.GeneID.quantile(threshold)] \
                .index.tolist()
            return selected_values_as_features

        def build_dataset(source_data, feature_type):

            settings = {'pathways': ('PathwayName', 0.9),
                        'chemicals': ('InteractionActions', 0.5)}
            grouping_col = settings[feature_type][0]

            features = select_features(source_data, grouping_col, settings[feature_type][1])
            selected_features_df = source_data[source_data[grouping_col].isin(features)]
            other_features = source_data[-source_data[grouping_col].isin(features)]

            if feature_type == 'pathways':
                grouped_features_df = selected_features_df.groupby(['GeneSymbol', grouping_col]).GeneID.nunique()
                other_features = other_features.groupby('GeneSymbol').GeneID.nunique()
                d = {'GeneSymbol': "Gene_Name", 'GeneID': 'otherPathways'}
                values_name = "GeneID"
            else:
                grouped_features_df = selected_features_df.groupby(['GeneSymbol', grouping_col]).ChemicalName.count()
                other_features = other_features.groupby('GeneSymbol').ChemicalName.count()
                unique_interactions = source_data.groupby('GeneSymbol').InteractionActions.nunique()
                other_features = pd.merge(unique_interactions, other_features, how='outer', on='GeneSymbol')
                d = {'ChemicalName': "otherInteractionsCount", 'GeneSymbol': "Gene_Name",
                     "InteractionActions": 'uniqueInteractions'}
                values_name = "ChemicalName"
            grouped_features_df = grouped_features_df.to_frame()
            grouped_features_df.reset_index(inplace=True)

            grouped_features_df = grouped_features_df.pivot_table(index='GeneSymbol', values=values_name,
                                                                  columns=grouping_col)

            total_features = pd.merge(grouped_features_df, other_features, how='outer', on='GeneSymbol')
            total_features.reset_index(inplace=True)
            total_features.rename(columns=d, inplace=True)
            total_features.columns = ['ctd_' + x if x != 'Gene_Name' else x for x in total_features.columns]
            # print(total_features.head())

            return total_features

        print("\n>> Compiling CTD chemical-genes features...")
        gene_df = pd.read_csv(self.cfg.data_dir / 'PHAROS/pharos_GF_wINDEX.csv', index_col=None)
        gene_df = gene_df[['Gene_Name']]
        gene_df.drop_duplicates(inplace=True)

        pathways = pd.read_csv(self.cfg.data_dir / 'CTD/processed_CTD_genes_pathways.csv.gz', index_col=False)
        chemicals = pd.read_csv(self.cfg.data_dir / 'CTD/processed_CTD_chem_gene_ixns.csv.gz', index_col=False)

        print("pathways")
        print(pathways.head())
        print("chemicals")
        print(chemicals.head())

        paths_df = build_dataset(pathways, 'pathways')
        chem_df = build_dataset(chemicals, 'chemicals')

        final_df = pd.merge(chem_df, paths_df, how='outer', on='Gene_Name')

        final_ctd_df = pd.merge(gene_df, final_df, how='left', on='Gene_Name')
        # print(final_ctd_df.head())
        print(final_ctd_df.shape)

        final_ctd_df.fillna(0, inplace=True)

        return final_ctd_df

    def process_string(self, seed_genes):
        '''
        Process features related to experimentally derived PPIs
        between genes being level-1 and level-2 neighbours
        :return: inweb_df
        '''

        print("\n>> Compiling STRINGdb features...")
        seed_genes = set(seed_genes)

        protein_df = self.compile_interaction_table('protein_links_interactions.tsv.gz', 'string_db', 'protein',
                                                    seed_genes)
        physical_df = self.compile_interaction_table('physical_links_interactions.tsv.gz', 'string_db', 'physical',
                                                     seed_genes)
        string_df = pd.merge(protein_df, physical_df, how='outer', left_on='Gene_Name', right_on='Gene_Name')

        datasets = ['protein', 'physical']

        for df in datasets:
            links_file = "string_db/processed_" + df + "_links.csv.gz"
            interactions_file = "string_db/" + df + "_links_interactions.tsv.gz"
            links = pd.read_csv(self.cfg.data_dir / links_file, index_col=None)
            interactions = pd.read_csv(self.cfg.data_dir / interactions_file, sep='\t', index_col=None)
            links.set_index('Gene_Name', inplace=True)
            interactions.set_index('Gene_Name', inplace=True)

            interactions['interacting_genes_score'] = links.groupby('Gene_Name')['combined_score'].apply(list)
            interactions['interacting_genes_score'] = interactions['interacting_genes_score'].apply(sum)

            interactions['overlap_genes_score'] = links[links.gene2.isin(seed_genes)].groupby('Gene_Name')[
                'combined_score'].apply(list)
            interactions['overlap_genes_score'].fillna("", inplace=True)
            interactions['overlap_genes_score'] = interactions['overlap_genes_score'].apply(sum)

            interactions["string_db_L1_" + df + "_seed_genes_overlap_weighted_score"] = interactions[
                                                                                            'overlap_genes_score'] / \
                                                                                        interactions[
                                                                                            'interacting_genes_score']

            interactions["string_db_L2_" + df + "_sseed_genes_overlap_weighted_score"] = interactions[
                'interacting_genes']. \
                apply(lambda x: interactions['overlap_genes_score'].loc[eval(x)].sum() /
                                interactions['interacting_genes_score'].loc[eval(x)].sum())

            interactions['interacting_genes_hmean_score'] = links.groupby('Gene_Name')['combined_score'].apply(list)
            interactions['interacting_genes_hmean_score'] = interactions['interacting_genes_hmean_score'].apply(hmean)
            interactions['overlap_genes_hmean_score'] = links[links.gene2.isin(seed_genes)].groupby('Gene_Name')[
                'combined_score'].apply(list)
            interactions['overlap_genes_hmean_score'].loc[interactions['overlap_genes_hmean_score'].isnull()] = \
                interactions['overlap_genes_hmean_score'].loc[interactions['overlap_genes_hmean_score'].isnull()].apply(
                    lambda x: [1])
            interactions['overlap_genes_hmean_score'] = interactions['overlap_genes_hmean_score'].apply(hmean)
            interactions["string_db_L1_" + df + "_sseed_genes_overlap_hmean_score"] = interactions[
                                                                                          'overlap_genes_hmean_score'] / \
                                                                                      interactions[
                                                                                          'interacting_genes_hmean_score']

            #         interactions["String_L2_seed_genes_overlap_hmean_score"] = interactions['interacting_genes'].\
            #             apply(lambda x: hmean(list(links[links.gene2.isin(seed_genes) & links.loc[eval(x)]].combined_score))/hmean(list(links.loc[eval(x)].combined_score)))

            interactions_df = interactions[["string_db_L1_" + df + "_seed_genes_overlap_weighted_score",
                                            "string_db_L2_" + df + "_sseed_genes_overlap_weighted_score", \
                                            "string_db_L1_" + df + "_sseed_genes_overlap_hmean_score"]]
            #         , "string_db_L2_seed_genes_overlap_hmean_score"]]

            string_df = pd.merge(string_df, interactions_df, how='outer', left_on='Gene_Name', right_on='Gene_Name')

        string_df.fillna(0, inplace=True)

        return string_df

    def process_omim_disease_count(self):
        '''
        Process data with unique number of diseases per gene derived computed based on OMIM data
        :return: disease_count_df
        '''
        print("\n>> Compiling OMIM unique disease counter...")

        # disease_count_df = pd.read_csv(self.cfg.data_dir / 'omim/ndiseases_perc75.txt', index_col=0)
        # # disease_count_df = pd.read_csv(self.cfg.data_dir / 'omim/ndiseases_perc50.txt', index_col=0)
        disease_count_df = pd.read_csv(self.cfg.data_dir / 'omim/ndiseases_pvalues.txt', index_col=0)

        disease_count_df = disease_count_df.rename(columns={'gene': "Gene_Name", "ndiseases": "OMIM_uniq_diseases"})
        disease_count_df.drop_duplicates(inplace=True)

        return disease_count_df

    ############ NEW
    def run_all(self, labeler):

        if len(set(['pharos', 'inter_pro']) & set(self.cfg.data_source)) > 0:

            druggability_datasets = []

            ### Pharos data processing
            if 'pharos' in self.cfg.data_source:

                pharos_datasets = []

                pharos_df = self.process_pharos()
                pharos_datasets.append(pharos_df)
                print('PHAROS:', pharos_df.shape)

                dgidb_df = self.process_dgidb()
                pharos_datasets.append(dgidb_df)
                print('DGIdb:', dgidb_df.shape)

                seed_genes = labeler.provide_seed_genes()
                reactome_df = self.process_reactome(seed_genes)
                pharos_datasets.append(reactome_df)

                inweb_df = self.process_inweb(seed_genes)
                pharos_datasets.append(inweb_df)
                print('InWeb_IM:', inweb_df.shape)

                ctd_df = self.process_ctd()
                pharos_datasets.append(ctd_df)
                print('CTD:', ctd_df.shape)

                if not hasattr(self.cfg, 'debug_mode'):
                    string_df = self.process_string(seed_genes)
                    pharos_datasets.append(string_df)
                    print('STRINGdb:', string_df.shape)
                else:
                    pass

                disease_count_df = self.process_omim_disease_count()
                pharos_datasets.append(disease_count_df)
                print('OMIM unique disease count:', disease_count_df.shape)

                print("\n>> Merging all data frames together...")
                druggability_features_df = reduce(lambda df1, df2: pd.merge(df1, df2, how='left', on='Gene_Name'), \
                                                  pharos_datasets)
                print('druggability features df shape: ', druggability_features_df.shape)

                druggability_features_df.drop_duplicates(subset='Gene_Name', inplace=True)
                druggability_datasets.append(druggability_features_df)

            ### InterPro data processing
            if hasattr(self.cfg, 'inter_pro'):
                interpro_features = []
                if 'domains' in self.cfg.inter_pro:
                    domain_df = self.process_interpro_domain()
                    interpro_features.append(domain_df)
                    print('InterPro Domains:', domain_df.shape)

                if 'families' in self.cfg.inter_pro:
                    family_df = self.process_interpro_family()
                    interpro_features.append(family_df)
                    print('InterPro Families:', family_df.shape)

                if 'super_families' in self.cfg.inter_pro:
                    super_family_df = self.process_interpro_super_family()
                    interpro_features.append(super_family_df)
                    print('InterPro Super Families:', super_family_df.shape)

                interpro_features_df = reduce(lambda df1, df2: pd.merge(df1, df2, how='left', on='Gene_Name'),
                                              interpro_features)
                interpro_features_df.drop_duplicates(subset=['Gene_Name'], inplace=True)
                druggability_datasets.append(interpro_features_df)
                print('InterPro shape: ', interpro_features_df.shape)

            druggability_df = reduce(lambda df1, df2: pd.merge(df1, df2, how='outer', on='Gene_Name'),
                                     druggability_datasets)
            print('final druggability df shape: ', druggability_df.shape)

            druggability_df.drop_duplicates(inplace=True)
            druggability_df.to_csv(self.cfg.druggability_feature_table, sep='\t', index=None)
            print("Saved to {0}".format(self.cfg.druggability_feature_table))

            print(druggability_df.shape)

        else:
            print("No druggability related data sources included")


if __name__ == '__main__':
    out_dir = sys.argv[1]
    cfg = Config(out_dir)
    cfg.tier_tag = ['Tier 1']
    proc = ProcessDruggabilityFeatures(cfg)
    proc.run_all()
