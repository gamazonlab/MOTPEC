# MOTPEC
## a Multi-Omics-based framework bridges Transcriptome between PEripheral and Central tissues

MOTPEC trains a prediction model using lasso regression based on peripheral blood profiles to predict the gene expression of central tissues. The peripheral blood profiles include expression, eSNP, splicing, APA event, co-expression module and transcriptional factors. For each gene in each tissue, MOTPEC predicts its tissue-specific expression.The pearson's r between observed expression and predicted expression is employed to evaluate the accuracy of model. Integrating TWAS, compare DGE of observed expression and predicted expression to estimate the utility of prediction model.

## data availability
The data used in this study is freely available for download in the GTEx portal, https://www.gtexportal.org/

## workdir
workdir: '~/MOTPEC'
raw_data: '~/MOTPEC/data/raw_data'
input_data: '~/MOTPEC/data/input'
output_data: '~/MOTPEC/data/output'

## main script
Codes of MOTPEC are as followes:

###get all input file
input files include expression, genetics variants, splicing, APA event, transcriptional factors list, demography variables, BMI phenotype. All these files should be kept in input. 
1) expression: normal expression of 48 tissues should be kept in '~/MOTPEC/data/raw_data/normal_exp', and formatted expression data will be store in '~/MOTPEC/data/iutput/exp_matrix.Rdata'
2) partial whole blood profiles include blood.pca, blood.tf.pca, modules_pca. transcriptional factors list should be kept in '~/MOTPEC/data/raw_data/TF.txt', and formatted data will be stored in '~/MOTPEC/data/iutput/WB.Rdata'
3) splicing: splicing file should be kept in '~/MOTPEC/data/raw_data/Whole_Blood.v8.leafcutter_phenotypes.bed.gz', and formatted data will be stored in '~/MOTPEC/data/iutput/splicing.Rdata'.
4) APA event: apa event file should be kept in '~/MOTPEC/data/raw_data/Whole_Blood_All_PDUIs.zip', samples attributes should be kept in '~/MOTPEC/data/raw_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', and formatted data will be stored in '~/MOTPEC/data/iutput/apa_pc_data.Rdata'.
5) demography variables: demography variables 'sample_ga.txt' should be kept in '~/MOTPEC/data/raw_data/sample_ga.txt', and formatted data will be stored in '~/MOTPEC/data/iutput/sam_demo.Rdata'
6) BMI phenotype: bmi phenotype 'GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt' should be kept in 'MOTPEC/data/raw_data/', and formatted data will be stored in '~/MOTPEC/data/iutput/bmi.Rdata'.
7) other file: loci file 'gencode.v32.GRCh37.txt' should be kept in '~/MOTPEC/data/iutput/'
The script 'get_input_data.R' is used to get all input files.

###construct prediction model
MOTPEC uses Lasso regression model under a 5-fold cross-validation to train prediction model, get prediction and evaluate accuracy. The results will be saved in '~/MOTPEC/data/output/prediction_model_rs/' by tissue, including baseline's accuracy, MOTPEC's accuracy, MOTPEC's prediction, MOTPEC's beta. The script 'get_input_file.R' is used to realize above work, and 'functions.R' saves some functions used in 'get_input_file.R'

###format predicted results
The script 'format_prediction.R' will format the prediction model results, and 'pcc_2models', 'exp_pred', 'all_model_beta', 'sam_gene_time' will be stored in '~/MOTPEC/data/output/'

###do TWAS
PrediXcan is employed to do TWAS in this study. The code is in folder 'do_TWAS'. The required data should be kept in '~/MOTPEC/data/raw_data/' beforehand. The TWAS results will be stored in '~/MOTPEC/data/output/twas_rs'.

###do DGE
A linear regression equation after confonding factors adjusted is used to estimate the association between expression and BMI. MOTPEC does this on observed expression, predicted expression and baseline expression respectively. The script 'do_DGE.R' is used to realize above work, and DGE results will be stored in '~/MOTPEC/data/output/b_beta.Rdata'

## Required Packages
For the code the following R-packages needs to be installed 
'data.table', 'WGCNA', 'readr', 'ggplot2', 'glmnet', 'caret', 'foreach', 'doMC'
