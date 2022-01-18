# -*- coding: utf-8 -*-

from supervenn import supervenn
import matplotlib.pyplot as plt
import ast
import sys

def read_data(model):
    therapy_areas_sm = {}
    therapy_areas_ab = {}
    therapy_areas_other = {}
    
    fin = open('./ot_clinical_data_'+model+'.tsv','r')
    
    line = fin.readline()
    
    for line in fin:
        row = line.split('\t')
        tas = ast.literal_eval(row[7])
        
        if row[4] == 'Small molecule':
            for ta in tas:
                if ta not in ['measurement', 'biological process', 'phenotype']:
                    if ta in therapy_areas_sm.keys():
                        therapy_areas_sm[ta].add(row[1])
                    else:
                        therapy_areas_sm[ta] = set([row[1]])
        elif row[4] == 'Antibody':
            for ta in tas:
                if ta not in ['measurement', 'biological process', 'phenotype']:
                    if ta in therapy_areas_ab.keys():
                        therapy_areas_ab[ta].add(row[1])
                    else:
                        therapy_areas_ab[ta] = set([row[1]])
        else:
            for ta in tas:
                if ta not in ['measurement', 'biological process', 'phenotype']:
                    if ta in therapy_areas_other.keys():
                        therapy_areas_other[ta].add(row[1])
                    else:
                        therapy_areas_other[ta] = set([row[1]])
                
    fin.close()

    return therapy_areas_sm, therapy_areas_ab, therapy_areas_other


def plot_therapy_areas(therapy_areas_sm, therapy_areas_ab, therapy_areas_other):
    plt.figure(figsize=(16, 8))
    supervenn(list(therapy_areas_sm.values()),list(therapy_areas_sm.keys()),chunks_ordering='minimize gaps',sets_ordering='size',min_width_for_annotation=4)
    
    plt.grid(b=None)
    plt.xlabel('Top DrugnomeAI Genes')
    plt.ylabel('Therapeutic Areas')
    plt.title('Small Molecule (DrugnomeAI-'+model+')', pad=100)
    
    plt.tight_layout()
    plt.savefig('./therapy_area_sm_'+model+'.png')
    plt.savefig('./therapy_area_sm_'+model+'.pdf')
    
    plt.figure(figsize=(16, 8))
    supervenn(list(therapy_areas_ab.values()),list(therapy_areas_ab.keys()),chunks_ordering='minimize gaps',sets_ordering='size',min_width_for_annotation=4)
    
    plt.grid(b=None)
    plt.xlabel('Top DrugnomeAI Genes')
    plt.ylabel('Therapeutic Areas')
    plt.title('Monoclonal Antibody (DrugnomeAI-'+model+')', pad=100)
    
    plt.tight_layout()
    plt.savefig('.therapy_area_ab_'+model+'.png')
    plt.savefig('.therapy_area_ab_'+model+'.pdf')
    
    plt.figure(figsize=(16, 8))
    supervenn(list(therapy_areas_other.values()),list(therapy_areas_other.keys()),chunks_ordering='minimize gaps',sets_ordering='size',min_width_for_annotation=4)
    
    plt.grid(b=None)
    plt.xlabel('Top DrugnomeAI Genes')
    plt.ylabel('Therapeutic Areas')
    plt.title('Other Molecules (DrugnomeAI-'+model+')', pad=100)
    
    plt.tight_layout()
    plt.savefig('./therapy_area_other_'+model+'.png')
    plt.savefig('./therapy_area_other_'+model+'.pdf')


if __name__ == '__main__': 
    model = sys.argv[1]
    
    therapy_areas_sm, therapy_areas_ab, therapy_areas_other = read_data(model)
    plot_therapy_areas(therapy_areas_sm, therapy_areas_ab, therapy_areas_other)
    

    