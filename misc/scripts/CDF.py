import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})

def plot_cdf(models):
    interval_index = list(range(0,20))
    fig, ax = plt.subplots(figsize=(8,8))
    
    for model in models:
        interval = []
        count = []
        
        fin = open('./clinic_lower_genes_fisher_'+model+'.tsv','r')
        
        line = fin.readline()
        
        for line in fin:
            row = line.split('\t')
            interval.append(row[0])
            count.append(int(row[1]))
            
        
        count = np.array(count)
        pdf = count / sum(count)
        cdf = np.cumsum(pdf)
        
        ax.plot(interval_index, cdf, label=model+'-based scores')
        
        fin.close()
    
    ax.set_title('CDF of genes supported by clinical trials evidence\nper rank interval')
    ax.set_xticks(np.arange(0,len(interval)))
    ax.set_xticklabels(interval, rotation=90)
    ax.set_xlabel('Gene ranking intervals')
    ax.set_ylabel('Percentage of genes supported by\nclinical evidence')
    
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('./CDF.pdf')


if __name__ == '__main__': 
    models = ['Tclin','Tier1']
    
    plot_cdf(models)