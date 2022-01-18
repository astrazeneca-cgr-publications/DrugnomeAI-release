import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import random
import sys

from scipy.integrate import simpson

plt.rcParams.update({'font.size': 24})

class StepwiseHypergeom:

    def __init__(self, full_ranking_files, external_successes_file, bucket_size, bucket='', mol_type='', full_ranking_top_perc=-1, successes_top_perc=-1, max_x_lim=-1, max_y_lim=-1):
        
        """
            This class calculates the stepwise hypergeometric enrichment between a full ranked list of genes and a set of external files containing ranked 'successes'.
            
            - full_ranking_file: 
                full reference gene ranking; file should be provided as a single-column ranked list of genes (withouth column header)  
           
            - external_successes_files: 
                list of external files containing the 'successes'; each file should be provided as a single-column ranked list of genes (withouth column header)
            
            - full_ranking_top_perc: 
                top percentage of full_ranking_file's hits to restrict the enrichment analysis into (e.g 0.05); if -1, the full list is used
                
            - successes_top_perc:
                top percentage of each of external_successes_files' hits to restrict the enrichment analysis into (e.g 0.05); if -1, the full list is used
                
            - max_x_lim: 
                x-axis limit for the hypergeometric enrichment plots
            
            - max_y_lim: 
                y-axis limit for the hypergeometric enrichment plots
        """
        
        #self.full_ranking_df = self.read_full_ranking(full_ranking_file)
        self.full_ranking_files = full_ranking_files
        
        self.external_successes = pd.read_csv(external_successes_file, header=None)
        
        self.full_ranking_top_perc = full_ranking_top_perc
        self.successes_top_perc = successes_top_perc
                
        self.bucket_size = bucket_size
        self.bucket = bucket
        self.mol_type = mol_type
        
        # plotting parameters
        self.figure_name = 'both'
        self.max_x_lim = max_x_lim
        self.max_y_lim = max_y_lim
        
    
    def calc_phred_score(self, pval):
        phred = None
        
        # Set an upper cap of 500 for all Phred scores
        if pval == 0:
            phred = 500
        else:
            phred = -10 * np.log10(pval)
            if phred > 500:
                phred = 500
        
        return phred



    def read_full_ranking(self, full_ranking_file):

        full_ranking_df = pd.read_csv(full_ranking_file, header=None)
        full_ranking_df.columns = ['Gene_Name']
        full_ranking_df['full_ranking'] = range(1, full_ranking_df.shape[0]+1)

        return full_ranking_df



    def calc_hypergeom_pvals(self, full_ranking_df):

        external_successes = self.external_successes.copy()
        # external top predictions
        external_successes.columns = ['Gene_Name']
        #print("\n external ranking:", external_successes.head())
        #print(external_successes.shape)


        # Subset top external hits to then overlap with the full rankings df
        if self.successes_top_perc != -1:
            #external_successes = external_successes.iloc[:int(external_successes.shape[0] * self.successes_top_perc), :]
            external_successes = external_successes.iloc[:self.bucket_size, :]
        external_top_genes = external_successes['Gene_Name'].tolist()
        print('- Total number of Successes:', len(external_top_genes))




        # Restrict full ranking to top perc results if 'full_ranking_top_perc' is specified
        if self.full_ranking_top_perc != -1:
            #print(full_ranking_df.shape)
            #full_ranking_df = full_ranking_df.loc[:int(full_ranking_df.shape[0] * self.full_ranking_top_perc)]
            full_ranking_df = full_ranking_df.iloc[:992]
            #print(full_ranking_df.shape)
            print('- Restrict full ranking to top % results: ' + str(int(self.full_ranking_top_perc*100)) + '% (' + str(full_ranking_df.shape[0]) + ')' )
        

        M = full_ranking_df.shape[0]
        print('- Population Size:', M)
        #print('\n Full ranking:\n', full_ranking_df.head())




        full_ranking_df = full_ranking_df.loc[full_ranking_df['Gene_Name'].isin(external_top_genes)]
        full_ranking_df.reset_index(drop=True, inplace=True)
        #print('\n', full_ranking_df.head())
        #print('\n', full_ranking_df.tail())

        
        n = full_ranking_df.shape[0]
        print('\nTotal number of Successes in Full ranking (overlap):', n)





        # ************* Hypergeometric Test *************
        hypergeom_pvals = []
        hypergeom_ordered_genes = []

        for x in range(full_ranking_df.shape[0]):

            #print(full_ranking_df.iloc[x, 0])

            N = full_ranking_df.iloc[x, full_ranking_df.shape[1]-1]
            #print(x, N)


            cur_pval = hypergeom.sf(x - 1, M, n, N)
            #print(cur_pval)

            hypergeom_pvals = hypergeom_pvals + [cur_pval]

            cur_gene = full_ranking_df.loc[x, 'Gene_Name']
            hypergeom_ordered_genes = hypergeom_ordered_genes + [cur_gene]
        # ***********************************************


        #min_pval = min(hypergeom_pvals)
        hypergeom_pvals = [self.calc_phred_score(pval) for pval in hypergeom_pvals]
        print('\nMax phred:', max(hypergeom_pvals))
        print('Avg phred:', np.mean(hypergeom_pvals))

        return hypergeom_pvals




    def plot_hypergeom_stepwise(self, hypergeom_pvals, ax, label, color='#33a02c', signif_thres=0.05, plot_pval_thres=True):

        if plot_pval_thres:
            signif_thres = self.calc_phred_score(signif_thres)
            ax.axhline(y=signif_thres, linestyle='--', color='red', label='p-val: 0.05')


        ax.plot(hypergeom_pvals, color=color,
                label=label,
                linewidth=1)

        ax.set_xlim(left=-0.5)

        y_label = 'Phred score from hypergeometric test\n(significance increasing in positive direction)'
        ax.set_ylabel(y_label, fontsize=14)

        ax.legend(bbox_to_anchor=(1.32, 1), fontsize=12, loc='upper right', framealpha =0.6)
        for line in ax.legend().get_lines():
            line.set_linewidth(2.0)


    def get_area_under_curve(self,yvals):
        # - Calc. AUC based on Simpson's area / trapezoid area
        # - Restrict area above the line: y=phred(0.05)
        
        phred = self.calc_phred_score(0.05)
        yvals = np.array(yvals)
        
        diff = yvals - phred
        curve_values = np.maximum(diff, 0) #negative values are replaced with 0
        curve_auc = simpson(curve_values, dx=1)
        
        return curve_auc
    
        
    def main(self):
        print("hello")
        hypergeom_results = {}
            
        # Calculate hypergeometric enrichments for each external successes file
        for f in self.full_ranking_files:
            #print(f)
            full_ranking_df = self.read_full_ranking(f)
            ranking_file = f.split('/')[-1]
            print('\n\n\n> ', ranking_file + ':', full_ranking_df.shape)    

            hypergeom_results[ranking_file] = self.calc_hypergeom_pvals(full_ranking_df)


        # # Calculate hypergeometric enrichments against a random gene ranking
        random_full_ranking = full_ranking_df.sample(frac=1, random_state=0)
        random_full_ranking.reset_index(drop=True, inplace=True)
        random_full_ranking['full_ranking'] = range(1, full_ranking_df.shape[0]+1)
        #print(random_full_ranking.head())
        print('\n\n\n> Random:', random_full_ranking.shape)   
        random_hypergeom_pvals = self.calc_hypergeom_pvals(random_full_ranking)
        
        print(sorted(hypergeom_results, key=lambda k: len(hypergeom_results[k]), reverse=True))        


        # ------------------------- Start plotting ----------------------------
        print('\n> Plotting hypergeom. enrichment plots...')
        fig, ax = plt.subplots(figsize=(12, 8))

        colors = ['#2171b5', '#238b45', '#fe9929', '#fff7bc']
        cnt = 0
        for external_dataset in sorted(hypergeom_results, key=lambda k: len(hypergeom_results[k]), reverse=True):  #hypergeom_results.items():
        
            hypergeom_pvals = hypergeom_results[external_dataset]
            
            #calculate area under curve
            auc = self.get_area_under_curve(hypergeom_pvals)
            print(external_dataset+" auc= ",auc)
            
            self.plot_hypergeom_stepwise(hypergeom_pvals, ax, label=external_dataset.split('.')[0], color=colors[cnt], plot_pval_thres=False)
            cnt += 1

        self.plot_hypergeom_stepwise(random_hypergeom_pvals, ax, label='random', color='black')
        
        # Set axis limits
        if self.max_x_lim != -1:
            ax.set_xlim(0, self.max_x_lim)
        else:
            ax.set_xlim(0, max(len(l) for k, l in hypergeom_results.items()))
        if self.max_y_lim != -1:
            ax.set_ylim(0, self.max_y_lim)

        # Define filename top perc. identifiers
        if self.successes_top_perc == -1:
            success_top_str = 'all'
        else:
            success_top_str = 'top' + str(int(self.successes_top_perc * 100)) + 'perc'
                 
        if self.full_ranking_top_perc == -1:
            full_top_str = 'all'
        else:
            full_top_str = 'top' + str(int(self.full_ranking_top_perc * 100)) + 'perc'
        
        if self.mol_type == 'Kings':
            plt.title('Kings et al')
        else:
            plt.title(self.mol_type + ' ' + self.bucket)
            
        # Save fig to pdf
        fig.savefig(self.mol_type+'_External_' + self.bucket + '-vs-Full_ranking_' + self.figure_name + '_' + full_top_str + '.pdf', bbox_inches='tight')
        fig.savefig(self.mol_type+'_External_' + self.bucket + '-vs-Full_ranking_' + self.figure_name + '_' + full_top_str + '.png', bbox_inches='tight')
        plt.close()
        
        print('\n> Done.')