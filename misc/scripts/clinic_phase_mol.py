import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import sys


def read_data(model):
    gene_list = pd.read_csv('./GradientBoostingClassifier.drugnome_ai_predictions_'+model+'.csv', sep=",")
    print(gene_list.head())
    
    glist = gene_list['Gene_Name'].squeeze().str.strip().tolist()
    top_genes = glist[0:int(len(glist)*0.05)]
    print(top_genes[:10])
    
    ot_gene_list = pd.read_csv('./ot_clinical_data_'+model+'.tsv', sep="\t")
    print(ot_gene_list.head())
    
    return top_genes, ot_gene_list


def analyze_phase_molecule_type(top_genes, ot_gene_list):
    ot_list = ot_gene_list['target'].unique()
    
    no_clinic_genes = len(top_genes) - len(ot_list)
    print(no_clinic_genes)
    
    clinic = {}
    clinic['sm'] = [0 for x in range(5)]
    clinic['ab'] = [0 for x in range(5)]
    clinic['other'] = [0 for x in range(5)]
    
    sm_genes = set()
    
    for x in [4,3,2,1,0]:
        sm_list = ot_gene_list.loc[(ot_gene_list['mol_type'] == 'Small molecule') & (ot_gene_list['phase'] == x)]
        sm_set = set(sm_list['target'].unique()) - sm_genes
        clinic['sm'][x] = len(sm_set)
        sm_genes = sm_genes | sm_set
        
    ab_genes = set()
    
    for x in [4,3,2,1,0]:
        ab_list = ot_gene_list.loc[(ot_gene_list['mol_type'] == 'Antibody') & (ot_gene_list['phase'] == x)]
        ab_set = set(ab_list['target'].unique()) - ab_genes
        clinic['ab'][x] = len(ab_set)
        ab_genes = ab_genes | ab_set
    
    other_genes = set()
    
    for x in [4,3,2,1,0]:
        other_list = ot_gene_list.loc[(~ot_gene_list['mol_type'].isin(['Small molecule', 'Antibody'])) & (ot_gene_list['phase'] == x)]
        other_set = set(other_list['target'].unique()) - other_genes
        clinic['other'][x] = len(other_set)
        other_genes = other_genes | other_set

    return clinic, no_clinic_genes


def draw_bargraph(model, clinic, no_clinic_genes):
    plt.rcParams.update({'font.size': 18})

    clinic['no_clinic'] = [0,0,0,0,no_clinic_genes]
    labels = ['Small Molecule','Antibody','Other','No clinical data']
    x = np.arange(len(labels))
    width = 0.2

    fig, ax = plt.subplots()
    fig.set_figheight(8)
    fig.set_figwidth(8)
    
    ax.bar(x - (width*2), [clinic['sm'][4],clinic['ab'][4],clinic['other'][4],clinic['no_clinic'][0]], width, label='Phase IV')
    ax.bar(x - width, [clinic['sm'][3],clinic['ab'][3],clinic['other'][3],clinic['no_clinic'][0]], width, label='Phase III')
    ax.bar(x, [clinic['sm'][2],clinic['ab'][2],clinic['other'][2],clinic['no_clinic'][0]], width, label='Phase II')
    ax.bar(x + width, [clinic['sm'][1],clinic['ab'][1],clinic['other'][1],clinic['no_clinic'][0]], width, label='Phase I')
    ax.bar(x + width, [clinic['sm'][0],clinic['ab'][0],clinic['other'][1],clinic['no_clinic'][4]], width, label='No Clinical Data')
    
    ax.set_ylabel('Number of genes')
    ax.set_title('Number of top ranked genes per molecule type\nand clinical trial phase (DrugnomeAI-'+model+')')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    plt.xticks(rotation=45, ha='right')
    ax.tick_params(axis=u'both', which=u'both',length=0)
    
    ax.legend()
    fig.tight_layout()
    
    plt.show()
    plt.savefig('./phase_mol_bar_'+model+'.pdf')


if __name__ == '__main__': 
    model = sys.argv[1]
    
    top_genes, ot_gene_list = read_data(model)
    clinic, no_clinic_genes = analyze_phase_molecule_type(top_genes, ot_gene_list)
    draw_bargraph(model, clinic, no_clinic_genes)
    
    