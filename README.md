-- Development of machine learning model integrating transcriptional and biological interaction data to predict non-small cell lung cancer using LASSO. --

* This repository contains all information to develope a LASSO model for
predicting Non-small cell lung cancer (NSCLC). *

It contains following steps:
1. Identification of Differentially Expressed Genes (DEGs)
    - Lung specific gene expression data was selected TCGA-TARGET-GTEx.
    - Ma'ayan lab's Appyter bulk RNA-seq analysis pipeline was used for DEGs.
        Jupyter Notebook "Xena_DE_Analysis_Pipeline.ipynb" for DEGs is provided.

2. Identification of biologically important nodes in the network
    - 40 genes identified by analyzing the DEGs and interactome.

3. LASSO model development 
    - Lung specific Gene Expression (RSEM norm_count) data was downloaded from UCSC Xena (https://xena.ucsc.edu/). 
    - The dataset consists of expression values of 40 genes identified in step (2), from 1013 samples of lung cancer and 397 samples as control, was selected.
    - This dataset was divided into 80% training dataset (TR contains 1128 samples) and 20% as independent test dataset 1 (TD1 contains 282 samples).
    - The LASSO model was developed using 10-fold cross-validations(cv) on the TR dataset, and performance was checked on test dataset TD1
    - The performance of LASSO model was checked on three independent datasets:
    GSE18842 (Sanchez-Palencia, Gomez-Morales et al. 2011), GSE27262 (Wei, Juan et al. 2012), and GSE19804 (Lu, Tsai et al. 2010).
    - Programe file (LASSO_Classification_NSCLC.R) is provided.
    - Input file 1: "GeneExp.txt", (Gene Expression data downloaded from Xena).
    - Input file 2: "DEGs_degree.txt", List of selecetd genes (40) identified by DEGs and Network analysis.
    - Input file 3: ""Sample_Name_3.csv", sample name and sample type (Cancer vs Normal). 
