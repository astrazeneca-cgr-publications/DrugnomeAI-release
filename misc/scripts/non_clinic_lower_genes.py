# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import json
import sys

plt.rcParams.update({'font.size': 13})

def calculate_fisher(success_interval, genes_clinic_all, interval_size, size):
    a = success_interval
    b = interval_size - a
    c = genes_clinic_all - a
    d = size - interval_size - c

    oddsratio, pvalue = stats.fisher_exact([[a, c], [b, d]],alternative='two-sided')

    return oddsratio, pvalue


def plot_oddsratio(values, keys, figure_name, ev_type, log=False):
    y = np.arange(0,len(keys))
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.plot(values, y, 'bo')
    ax.hlines(y, 0, values, lw=2)
    
    ax.set_yticks(y)
    ax.set_yticklabels(keys)
    
    if log:
        ax.plot([0 for x in range(len(y))], y, 'r-.', label='odds ratio = 1')
        ax.set_xlabel('log(Odds Ratio)')
    else:
        ax.plot([1 for x in range(len(y))], y, 'r-.', label='odds ratio = 1')
        ax.set_xlabel('Odds Ratio')
    
    ax.set_ylabel('ranking intervals')
    plt.gca().invert_yaxis()
    plt.title('Enrichment of genes supported by ' + ev_type+'\nevidence per rank interval (DrugnomeAI-'+model+')')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('./'+figure_name+'.png')
    plt.savefig('./'+figure_name+'.pdf')
    

def plot_pvalue(values, keys, figure_name, ev_type, log=False):
    y = np.arange(0,len(keys))
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.plot(values, y, 'bo')
    ax.hlines(y, 0, values, lw=2)
    
    ax.set_yticks(y)
    ax.set_yticklabels(keys)
    
    if log:
        ax.plot([np.log(0.05) for x in range(len(y))], y, 'r-.', label='p-value = 0.05')
        ax.set_xlabel('log(p-value)')
    else:
        ax.plot([0.05 for x in range(len(y))], y, 'r-.', label='p-value = 0.05')
        ax.set_xlabel('p-value')
    
    ax.set_ylabel('ranking intervals')
    plt.gca().invert_yaxis()
    plt.title('Enrichment of genes supported by ' + ev_type+'\nevidence per rank interval (DrugnomeAI-'+model+')')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('./'+figure_name+'.png')
    plt.savefig('./'+figure_name+'.pdf')
    
    
def write_results(successes, fisher_pvalue, log_pvalue, fisher_odds_ratio, log_oddsratio, interval_labels, file_name):
    fout = open('./lower_genes_fisher_'+file_name+'.tsv','w')
    fout.write("Interval\tSuccess\tFisher's p-value\tFisher's log(p-value)\tFisher's odds-ratio\t Fisher's log(odds-ratio)\n")
    
    for x in range(len(interval_labels)):
        fout.write(interval_labels[x]+'\t'+str(successes[x])+'\t'+str(fisher_pvalue[x])+'\t'+str(log_pvalue[x])+'\t'+str(fisher_odds_ratio[x])+'\t'+str(log_oddsratio[x])+'\n')
    
    fout.close()
    

def read_data(model):
    genes = pd.read_csv('./GradientBoostingClassifier.drugnome_ai_predictions_'+model+'.csv', sep=",")
    glist = genes['Gene_Name'].squeeze().str.strip().tolist()
    size = len(glist)
    
    ot_genes = pd.read_csv("./ot_clinical_data_all_genes.tsv", sep="\t")
    ot_genes_set = set(ot_genes['target'].squeeze().str.strip().tolist())
    
    return ot_genes_set, size, glist
    

def read_non_clinical_evidence(evidence_type, ot_genes_set):
    evidence_table = {}
    
    for ev_type in evidence_type:
        evidence_table[ev_type] = set()
    
    fin = open('../21.02_evidence_data.json','r')
    
    for line in fin:
        data = json.loads(line)
        target = data['target']['gene_info']['symbol']
        
        if data['type'] == 'genetic_association':
            evidence_table['Genetic'].add(target)
        elif data['type'] == 'affected_pathway':
            evidence_table['Pathway'].add(target)
        elif data['type'] == 'animal_model':
            evidence_table['Animal Models'].add(target)
        elif data['type'] == 'rna_expression':
            evidence_table['RNA Expression'].add(target)
        elif data['type'] == 'somatic_mutation':
            evidence_table['Somatic Mutation'].add(target)
        elif data['type'] == 'literature':
            evidence_table['Literature'].add(target)
            
    inter = np.arange(0.0, 1.05, 0.05)
    
    ranks = []
    non_clinic_genes = []
    
    for x in range(1,len(inter)):
        top_genes = set(glist[int(size*inter[x-1]):int(size*inter[x])])
        c_genes = len(top_genes & ot_genes_set)
        non_clinic_genes.append(len(top_genes) - c_genes)
        ranks.append("{0:.0%}".format(inter[x-1])+"-"+"{0:.0%}".format(inter[x]))
    
    return evidence_table, ranks, non_clinic_genes


def draw_bargraph(model, ranks, non_clinic_genes):
    width = 1
    
    fig, ax = plt.subplots(figsize=(8,6))
    rects = ax.barh(ranks, non_clinic_genes, width)
     
    ax.set_ylabel('Genes Ranks')
    ax.set_xlabel('Number of genes')
    ax.invert_yaxis()
    ax.bar_label(rects, padding=3)
    
    ax.set_title('Number of genes not supported by clinical trials\nevidence per rank interval (DrugnomeAI-'+model+')')
    plt.tight_layout()
    plt.show()
    plt.savefig('./lower_genes_no_clinic_'+model+'.png')
    plt.savefig('./lower_genes_no_clinic_'+model+'.pdf')
            

def plot_non_clinical_evidence(model, ev_type, ranks, non_clinic_genes, sample_size):
    width = 1
        
    fig, ax = plt.subplots(figsize=(8,6))
    rects = ax.barh(ranks, non_clinic_genes, width)
    
    ax.set_ylabel('Genes Ranks')
    ax.set_xlabel('Number of genes')
    ax.invert_yaxis()
    ax.bar_label(rects, padding=3)
    
    ax.set_title('Number of genes supported by '+ev_type+'\nevidence per rank interval (DrugnomeAI-'+model+')')
    plt.xlim(0,non_clinic_genes[0]+100)
    plt.tight_layout()
    plt.show()
    plt.savefig('./lower_genes_'+ev_type+'_'+model+'.png')
    plt.savefig('./lower_genes_'+ev_type+'_'+model+'.pdf')
    
    successes = []
    pvalue = []
    oddsratio = []
    non_genes_clinic_all = sum(non_clinic_genes)
    
    for y in interval_index:
        odd, p = calculate_fisher(non_clinic_genes[y], non_genes_clinic_all, sample_size[y], size)
        successes.append(non_clinic_genes[y])
        oddsratio.append(odd)
        pvalue.append(p)
        
    return successes, oddsratio, pvalue

            
if __name__ == '__main__': 
    model = sys.argv[1]
    evidence_type = ['Genetic', 'Pathway', 'Animal Models', 'RNA Expression', 'Somatic Mutation', 'Literature']

    ot_genes_set, size, glist = read_data(model)
    evidence_table, ranks, non_clinic_genes = read_non_clinical_evidence(evidence_type, ot_genes_set)    
    draw_bargraph(model, ranks, non_clinic_genes)

    inter = np.arange(0.0, 1.05, 0.05)
    interval_index = list(range(0,20))
    interval_labels = ["0%-5%", "5%-10%", "10%-15%", "15%-20%", "20%-25%", "25%-30%", "30%-35%", "35%-40%", "40%-45%", "45%-50%", "50%-55%", "55%-60%", "60%-65%", "65%-70%", "70%-75%", "75%-80%", "80%-85%", "85%-90%", "90%-95%", "95%-100%"]
    
    for ev_type in evidence_type:
        ranks = []
        non_clinic_genes = []
        sample_size = []
        
        for x in range(1,len(inter)):
            top_genes = set(glist[int(size*inter[x-1]):int(size*inter[x])])
            sample_size.append(len(top_genes))
            evidence_genes = top_genes & evidence_table[ev_type]
            non_clinic_genes.append(len(evidence_genes))
            ranks.append("{0:.0%}".format(inter[x-1])+"-"+"{0:.0%}".format(inter[x]))
        
        successes, oddsratio, pvalue = plot_non_clinical_evidence(model, ev_type, ranks, non_clinic_genes, sample_size)
            
        plot_oddsratio(oddsratio, interval_labels, 'oddsratio_'+ev_type+'_'+model, ev_type, log=False)
        log_oddsratio = np.log(oddsratio)
        plot_oddsratio(log_oddsratio, interval_labels, 'oddsration_log_'+ev_type+'_'+model, ev_type, log=True)
        
        plot_pvalue(pvalue, interval_labels, 'pvalue_'+ev_type+'_'+model, ev_type, log=False)
        log_pvalue = np.log(pvalue)
        plot_pvalue(log_pvalue, interval_labels, 'pvalue_log_'+ev_type+'_'+model, ev_type, log=True)
        
        write_results(successes, pvalue, log_pvalue, oddsratio, log_oddsratio, interval_labels, ev_type+'_'+model)


