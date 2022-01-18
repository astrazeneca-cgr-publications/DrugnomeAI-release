# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

plt.rcParams.update({'font.size': 13})


def calculate_fisher(success_interval, genes_clinic_all, interval_size, size):
    a = success_interval
    b = interval_size - a
    c = genes_clinic_all - a
    d = size - interval_size - c

    oddsratio, pvalue = stats.fisher_exact([[a, c], [b, d]],alternative='two-sided')

    return oddsratio, pvalue


def plot_oddsratio(values, keys, figure_name, model, log=False):
    y = np.arange(0,len(keys))
    
    fig = plt.figure(figsize=(8,6))
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
    plt.title('Enrichment of PROTACtable genes\nper rank interval (DrugnomeAI-'+model+')')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('./'+figure_name+'_'+model+'.png')
    plt.savefig('./'+figure_name+'_'+model+'.pdf')
    
    
def plot_pvalue(values, keys, figure_name, model, log=False):
    y = np.arange(0,len(keys))
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.plot(values, y, 'bo')
    ax.hlines(y, 0, values, lw=2)
    
    ax.set_yticks(y)
    ax.set_yticklabels(keys)
    
    if log:
        ax.plot([np.log(0.05) for x in range(len(y))], y, 'r-.', label='pvalue = 0.05')
        ax.set_xlabel('log(p-value)')
    else:
        ax.plot([0.05 for x in range(len(y))], y, 'r-.', label='pvalue = 0.05')
        ax.set_xlabel('p-value')
    
    ax.set_ylabel('ranking intervals')
    plt.gca().invert_yaxis()
    plt.title('Enrichment of PROTACtable genes\nper rank interval (DrugnomeAI-'+model+')')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('./'+figure_name+'_'+model+'.png')
    plt.savefig('./'+figure_name+'_'+model+'.pdf')
    
    
def write_results(successes,fisher_pvalue, log_pvalue, fisher_odds_ratio, log_oddsratio, interval_labels):
    fout = open('./clinic_lower_genes_fisher_'+model+'.tsv','w')
    fout.write("Interval\tSuccess\tFisher's p-value\tFisher's log(p-value)\tFisher's odds-ratio\t Fisher's log(odds-ratio)\n")
    
    for x in range(len(interval_labels)):
        fout.write(interval_labels[x]+'\t'+str(successes[x])+'\t'+str(fisher_pvalue[x])+'\t'+str(log_pvalue[x])+'\t'+str(fisher_odds_ratio[x])+'\t'+str(log_oddsratio[x])+'\n')
    
    fout.close()
    

def read_data(model):
    genes = pd.read_csv('./GradientBoostingClassifier.drugnome_ai_predictions_'+model+'.csv', sep=",")
    glist = genes['Gene_Name'].squeeze().str.strip().tolist()
    size = len(glist)
    
    fin = open('./protac_testing.txt','r')
    protac_genes = fin.readlines()
    protac_genes = [gene.strip() for gene in protac_genes]
    protac_genes_set = set(protac_genes)
    
    inter = np.arange(0.0, 1.05, 0.05)
    
    ranks = []
    clinic_genes = []
    sample_size = []

    for x in range(1,len(inter)):
        top_genes = set(glist[int(size*inter[x-1]):int(size*inter[x])])
        sample_size.append(len(top_genes))
        c_genes = len(top_genes & protac_genes_set)
        clinic_genes.append(c_genes)
        ranks.append("{0:.0%}".format(inter[x-1])+"-"+"{0:.0%}".format(inter[x]))

    return ranks, clinic_genes, sample_size, size


def draw_bargraph(model,ranks,clinic_genes):
    width = 1
    
    fig, ax = plt.subplots(figsize=(8,6))
    rects = ax.barh(ranks, clinic_genes, width)
     
    ax.set_ylabel('Genes Ranks')
    ax.set_xlabel('Number of genes')
    ax.invert_yaxis()
    ax.bar_label(rects, padding=3)
    
    ax.set_title('Number of PROTACtable genes\nper rank interval (DrugnomeAI-'+model+')')
    plt.xlim(0,clinic_genes[0]+50)
    plt.tight_layout()
    plt.show()
    plt.savefig('./lower_genes_'+model+'.png')
    plt.savefig('./lower_genes_'+model+'.pdf')


if __name__ == '__main__': 
    model = 'PROTACs'
    
    ranks, clinic_genes, sample_size, size = read_data(model)
    draw_bargraph(model,ranks,clinic_genes)
    
    interval_index = list(range(0,20))
    interval_labels = ["0%-5%", "5%-10%", "10%-15%", "15%-20%", "20%-25%", "25%-30%", "30%-35%", "35%-40%", "40%-45%", "45%-50%", "50%-55%", "55%-60%", "60%-65%", "65%-70%", "70%-75%", "75%-80%", "80%-85%", "85%-90%", "90%-95%", "95%-100%"]
    successes = []
    pvalue = []
    oddsratio = []
    genes_clinic_all = sum(clinic_genes)
    
    for inter in interval_index:
        odd, p = calculate_fisher(clinic_genes[inter], genes_clinic_all, sample_size[inter], size)
        successes.append(clinic_genes[inter])
        oddsratio.append(odd)
        pvalue.append(p)
        
    plot_oddsratio(oddsratio, interval_labels, 'oddsratio', model, log=False)
    log_oddsratio = np.log(oddsratio)
    plot_oddsratio(log_oddsratio, interval_labels, 'oddsration_log', model, log=True)
    
    plot_pvalue(pvalue, interval_labels, 'pvalue', model, log=False)
    log_pvalue= np.log(pvalue)
    plot_pvalue(log_pvalue, interval_labels, 'pvalue_log', model, log=True)
    
    write_results(successes, pvalue, log_pvalue, oddsratio, log_oddsratio, interval_labels)
    
    