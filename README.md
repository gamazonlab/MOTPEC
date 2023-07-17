# MOTPEC
## a Multi-Omics-based framework bridges Transcriptome between PEripheral and Central tissues

MOTPEC trains a prediction model using lasso regression based on peripheral blood profiles to predict the gene expression of central tissues. The peripheral blood profiles include expression, eSNP, splicing, APA event, co-expression module and transcriptional factors. For each gene in each tissue, MOTPEC predicts its tissue-specific expression.The pearson's r between observed expression and predicted expression is employed to evaluate the accuracy of model. Integrating TWAS, compare DGE of observed expression and predicted expression to estimate the utility of prediction model.

## data availability
The data used in this study is freely available for download in the GTEx portal, https://www.gtexportal.org/

## workdir
workdir: ~/MOTPEC
raw_data: ~/MOTPEC/data/raw_data
input_data: ~/MOTPEC/data/input
output_data: ~/MOTPEC/data/output

## main script
Codes of MOTPEC are as followes:

### step1: get all input file
input files include expression, genetics variants, splicing, APA event, transcriptional factors list, demography variables, BMI phenotype. All these files should be kept in input. 
1) expression: normal expression of 48 tissues should be kept in ~/MOTPEC/data/raw_data/normal_exp, and formatted expression data will be store in ~/MOTPEC/data/iutput/exp_matrix.Rdata 
2) partial whole blood profiles include blood.pca, blood.tf.pca, modules_pca. Transcriptional factors list has been kept in ~/MOTPEC/data/raw_data/annotation/TF.txt, colors and modules have been kept in ~/MOTPEC/data/raw_data/annotation/co_exp_module.Rdata, and formatted data will be stored in ~/MOTPEC/data/iutput/WB.Rdata
3) splicing: splicing file should be kept in ~/MOTPEC/data/raw_data/Whole_Blood.v8.leafcutter_phenotypes.bed.gz, and formatted data will be stored in ~/MOTPEC/data/iutput/splicing.Rdata
4) APA event: apa event file should be kept in ~/MOTPEC/data/raw_data/Whole_Blood_All_PDUIs.zip, samples attributes should be kept in ~/MOTPEC/data/raw_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt, and formatted data will be stored in ~/MOTPEC/data/iutput/apa_pc_data.Rdata.
5) demography variables: demography variables sample_ga.txt should be kept in ~/MOTPEC/data/raw_data/sample_ga.txt, and formatted data will be stored in ~/MOTPEC/data/iutput/sam_demo.Rdata
6) BMI phenotype: bmi phenotype GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt should be kept in MOTPEC/data/raw_data/, and formatted data will be stored in ~/MOTPEC/data/iutput/bmi.Rdata.
7) other file: loci file gencode.v32.GRCh37.txt should be kept in ~/MOTPEC/data/iutput/
The script get_input_data.R is used to get all input files.

### examples of input file
1) exp_matrix: a list of 49 tissues' expression(include Whole Blood), each element is a matrix
   > exp.matrix[['Whole_Blood']][1:5,1:5]
           ENSG00000227232 ENSG00000238009 ENSG00000233750 ENSG00000268903 ENSG00000269981
GTEX-111YS      -1.2250236      -0.6733178       0.3014908       0.4171043       0.8916543
GTEX-1122O      -0.8533908       0.1217088       0.8587798       0.3766841       0.3093189
GTEX-1128S       0.4293614      -0.6546895      -0.4375687       0.8163377       1.1788697
GTEX-113IC       0.8533908      -1.1788697      -1.5580346      -2.6142683      -2.6142683
GTEX-113JC      -0.7064943      -0.8373702      -0.7702820      -0.9312701      -0.8215617
2) WB: a list include


### step2: construct prediction model
MOTPEC uses Lasso regression model under a 5-fold cross-validation to train prediction model, get prediction and evaluate accuracy. The results will be saved in ~/MOTPEC/data/output/prediction_model_rs/ by tissue, including baseline accuracy, MOTPEC's accuracy, MOTPEC's prediction, MOTPEC's beta. The script get_input_file.R is used to realize above work, and functions.R saves some functions used in get_input_file.R

### step3: format predicted results
The script format_prediction.R will format the prediction model results, and pcc_2models, exp_pred, all_model_beta, sam_gene_time will be stored in ~/MOTPEC/data/output/

### step4: do TWAS
PrediXcan is employed to do TWAS in this study. The code is in folder do_TWAS. The required data should be kept in ~/MOTPEC/data/raw_data/ beforehand. The TWAS results will be stored in ~/MOTPEC/data/output/twas_rs.

### step5: do DGE
A linear regression equation after confonding factors adjusted is used to estimate the association between expression and BMI. MOTPEC does this on observed expression, predicted expression and baseline expression respectively. The script do_DGE.R is used to realize above work, and DGE results will be stored in ~/MOTPEC/data/output/b_beta.Rdata

### step6: predict the tissue-specific expression

## Required Packages
For the code the following R-packages needs to be installed 
data.table, WGCNA, readr, ggplot2, glmnet, caret, foreach, doMC
