# DrugnomeAI: gene druggability


- [Introduction](#introduction) 
- [Installation](#installation) 
- [Run](#run) 
  - [drugnomeai](#drugnomeai)



Introduction
============
`DrugnomeAI` is an adapatation of `mantis-ml` that provides both disease-agnostic and disease-specific gene druggability framework, implementing stochastic semi-supervised learning on top of `scikit-learn` and `keras`/`tensorflow`.  
`Drugnome` takes its name from a contraction of the druggable genome.

Installation
============
**Requirements:** `Python3` (tested with v3.6.7)

```
git clone https://github.com/astrazeneca-cgr-publications/DrugnomeAI-release.git
sh global_init.sh
python setup.py install
```

<br>

---

In either case, it is highly recommended to **create a new virtual environment** (e.g. with `conda`) before installing `DrugnomeAI`:
```
conda create -n drugnome python=3.6
conda config --append channels conda-forge   	# add conda-forge in the channels list
conda activate drugnome			                  # activate the newly created conda environment
```

---

<br>


You may now call the following scripts from the command line:
- **`drugnomeai`**: run DrugnomeAI druggability score based on a dummy config file (`.yaml`)

Run each command with `-h` to see all available options.


<br><br>



Run
===

### Disease-agnostic analysis
In this mode, the only things required to be specified are output directory and either type of labels or a file with custom seed genes which are as follows:
<br>
 `-o OUTPUT_DIR` or `--output-dir OUTPUT_DIR` and `-p clin` or `--pharos-tag clin`
 <br>
 `-o OUTPUT_DIR` or `--output-dir OUTPUT_DIR` and `-t 1` or `--tier-tag 1`
 <br>
 `-o OUTPUT_DIR` or `--output-dir OUTPUT_DIR` and `-k seed_genes_file` or `--known-genes-file seed_genes_file`

#### disease-agnostic run examples
```
drugnomeai -o out/Genes_analysis -n 4 -k custom_seed_genes.txt -r pre
drugnomeai -o out/Tier1_run -f -t 1
drugnomeai -o analysis_output -r boruta -p tchem 
```



<br>
<hr>
<br>

### Disease-specific analysis
You need to provide a config file (`.yaml`) containing information about the diseases/phenotypes of interest. 
<br>


##### Required field:
- `Disease/Phenotype terms`: **terms that characterise a phenotype or disease of interest** (*free text*)


##### Optional fields:
- `Additional associated terms`: terms used along with `Disease/Phenotype terms` to extract additional disease/phenotype-associated features (*free text*)
- `Diseases/Phenotypes to exclude`: terms to exclude from disease/phenotype characterisation and feature selection (*free text*)


<br>


**Config examples**:
```
# Epilepsy_config.yaml
Disease/Phenotype terms: epileptic, epilepsy, seizure
Additional associated terms: brain, nerve, nervous, neuronal, cerebellum, cerebral, hippocampus, hypothalamus
Diseases/Phenotypes to exclude: 
```
```
# CKD_config.yaml
Disease/Phenotype terms: renal, kidney, nephro, glomerul, distal tubule 
Additional associated terms: 
Diseases/Phenotypes to exclude: adrenal
```

Other example config files can be found under [example-input](example-input) or `drugnome_ai/conf`. 

<br>
<br>
 

### More about DrugnomeAI
 
#### Labels
- There are two types of labels available in DrugnomeAI - Pharos labels and the ones representing Tiers. Only one type can be be used at once, however it is possible to select multiple labels of the same type for a single run. 
- The available options are the following:
    - Pharos labels:
         - `tclin`: Tclin
         - `tchem`: Tchem
         - `tbio`: Tbio
         - `tdark`: Tdark

    - Tiers labels:
         - `1`: Tier 1
         - `2`: Tier 2
         - `3B`: Tier 3A
         - `3A`: Tier 3B

Label option may be specified in the following way: `-p tchem` or -t 3A. For multiple choices, use ' ' (space) as a separator, e.g. `-t 1 2 3A` etc.
In the disease-specific run, labels of a given choice are fetched automatically by DrugnomeAI for the genes being associated with the specified disease.

#### Supervised learning models
- `DrugnomeAI` runs 8 different supervised models by default: Extra Trees, Random Forest, SVC, Gradient Boosting, XGBoost, Deep Neural Net, Stacking Classifier and Naive Bayes. 
- It is also possible to run `DrugnomeAI` with the `-f / --fast` option, which will force DrugnomeAI to train only 4 classifiers: `Extra Trees`, `Random Forest`, `SVC` and `Gradient Boosting`.
- Additionally, the user may explicitly specify which supervised models to be used for training via the `-s` option. The available model options are coded as follows:
  - `et`: Extra Trees
  - `rf`: Random Forest
  - `gb`: Gradient Boosting
  - `xgb`: XGBoost
  - `svc`: Support Vector Classifier
  - `dnn`: Deep Neural Net
  - `stack`: Stacking classifier
  - `nb`: Naive Bayes

Multiple models may be delimited using a ` ` (space) as a separator, e.g. `-s et`, `-s et stack gb` etc. 

#### Data sources and features
`DrugnomeAI` allows to make a choice about features which can be selected from different sets which are the following:
    
   - `pharos`: denotes all druggability related features including: Pharos, DGIdb, Reactome, InWeb, CTD, STRINGdb
   - `inter`: InterPro features related to protein families, domains and super-families 
  
 By default, DrugnomeAI runs with all the features listed above. However, the feature set can be modified with `-d` option, e.g. `-d inter`.
 
 The disease-specific analysis, on top of the abovementioned, by default leverages features being specifically related to diseases and they are following:
 GTEx, Protein Atlas, features derived from collapsing analysis of Human Protein Atlas tissue & RNA expression and MsigDB GO.  
 
 The InterPro features can also be selected (choices: dom, fam, sup) by user with the option `-x dom fam` or `--inter-pro dom fam` - note that the options are delimited with space.
   
 There is also an option to extend generic DrugnomeAI features with `-m`  or `--mantis-ml` argument. 
   The option `--mantis-ml` denotes attributes derived from `mantis-ml` which includes the following:
    ExAC features, Essential Mouse Genes ones, GnomAD, Genic Intolerance Scores, GWAS & MGI Essential features.

  Another available option is `-l` or `--genic-intol` which allows you to include solely Genic Intolerance Scores out of all the mantis-ml derived features.
    



<br><br>


drugnomeai
===========
You need to provide an output directory and label choice or these along with a config file (`.yaml`). 
<br>
You may also:
- run specific modules of `DrugnomeAI` (option `-r`, default: all) - the available options are:
    - `all`: all modules are run, meaning end-to-end analysis shall be performed
    - `pre`: exploratory data analysis shall be run incl. data cleaning & visualizations
    - `boruta`: feature selection with Boruta algorithms on the already prepared data
    - `pu`: semi-supervised learning run
    - `post`: postprocessing of the results yielded by semi-supervised learning
    - `post_unsup`: dimensionality reduction & clustering run on the results
    - `debug`: initial exploratory analysis and 1 iteration of semi-supervised learning will be run without Inweb and StringDB features creation
- choose what tier of druggability you would like to train on (`-t` option; default: 1,2,3A,3B)
- choose to train on pharos based labels (`-p` option; default: clin,chem,bio,dark)
- define the number of threads to use (`-n` option; default value: 4).
- define the number of stochastic iterations (`-i` option; default value: 10)
- provide a file with custom seed genes (`-k` option; file should contain new-line separated HGNC names without header; bypasses HPO)

disease-agnostic run
```
drugnomeai -o [output_dir] [-n nthreads] [-i iterations] -k [custom_seed_genes.txt] [-r pre]
or 
drugnomeai -o [output_dir] [-n nthreads] [-i iterations] [-r pre] -t [druggability tier]
or
drugnomeai -o [output_dir] [-n nthreads] [-i iterations] [-r pre] -p [pharos druggability label]
```

disease-specific run
```
drugnomeai -c [config_file] -o [output_dir] [-n nthreads] [-i iterations] [-k custom_seed_genes.txt] [-r pre]
or 
drugnomeai -c [config_file] -o [output_dir] [-n nthreads] [-i iterations] [-r pre] [-p pharos druggability label]
```

#### Example
```
drugnomeai -c CKD_config.yaml -o ./CKD-run -f -p tclin
drugnomeai -o ./CKD-run -f -k kcd_seed_genes.txt
drugnomeai -c Epilepsy_config.yaml -o /tmp/Epilepsy-testing -n 20 -r boruta -p tclin
```


#### `drugnomeai` Output
`DrugnomeAI` predictions for all genes and across all classifiers can be found at **`[output_dir]/Gene-Predictions`**. 
<br>
The `AUC_performance_by_Classifier.pdf` file under the same dir contains information about the AUC performance per classifier and thus informs about the best performing classifier.

Output figures from all steps during the `mantis-ml` run (e.g. *Exploratory Data Analysis/EDA, supervised-learning, unsupervised-learning*) can be found under **`[output_dir]/Output-Figures`**.

<br>
