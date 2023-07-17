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
   ```exp.matrix[['Whole_Blood']][1:5,1:5]```  
               ENSG00000227232 ENSG00000238009 ENSG00000233750 ENSG00000268903 ENSG00000269981  
GTEX-111YS      -1.2250236      -0.6733178       0.3014908       0.4171043       0.8916543  
GTEX-1122O      -0.8533908       0.1217088       0.8587798       0.3766841       0.3093189  
GTEX-1128S       0.4293614      -0.6546895      -0.4375687       0.8163377       1.1788697  
GTEX-113IC       0.8533908      -1.1788697      -1.5580346      -2.6142683      -2.6142683  
GTEX-113JC      -0.7064943      -0.8373702      -0.7702820      -0.9312701      -0.8215617  
2) WB: a list include: blood.pca, blood.tf.pca, colors, modules_pca, modules  
   ```WB[['blood.pca']][1:5,1:5]```  
                     PC1       PC2        PC3        PC4        PC5  
GTEX-111YS  -68.56067  18.53914  -15.42782  -3.938970  20.411347  
GTEX-1122O -102.96515 -51.41632   29.38948 -40.438315 -13.367777  
GTEX-1128S   62.51361 -22.31329   27.41144 -65.461205  10.169594  
GTEX-113IC   16.41480 -83.65493 -104.79319  -5.910071  -1.902242  
GTEX-113JC   38.33161 -12.52743   60.43982  -9.061731   5.730263  
   ```WB[['blood.tf.pca']][1:5,1:5]```  
                     PC1           PC2        PC3        PC4       PC5  
GTEX-111YS  15.956726   8.756085788  -4.649753  -3.299114 -4.446283  
GTEX-1122O  25.759200  -5.709656529  11.203563  -9.575776  1.353987  
GTEX-1128S -15.255152  -7.904850477   1.738098 -14.940524 -5.970758  
GTEX-113IC   9.744112 -36.092799242 -23.263840   2.306993  1.778211  
GTEX-113JC  -9.658695  -0.004755886  12.989872  -2.667919 -2.465581  
   ```WB[['colors']][1:5]```  
   ENSG00000227232 ENSG00000238009 ENSG00000233750 ENSG00000268903 ENSG00000269981   
         "grey"          "grey"     "turquoise"     "turquoise"     "turquoise"   
   ```WB[['modules']][['grey']][1:5]```  
   ENSG00000227232 ENSG00000238009 ENSG00000241860 ENSG00000279928 ENSG00000279457   
              1               2               7               8               9  
   ```WB[['modules_pca']][['grey']][1:5,1:5]```  
                     PC1        PC2       PC3       PC4        PC5  
GTEX-111YS -27.653046   3.697308 -8.865022  2.724629 -19.317787  
GTEX-1122O -32.551785 -12.059564 32.652657  6.593311   1.591776  
GTEX-1128S  16.218895  -3.304558 21.695663 34.882897 -16.463806  
GTEX-113IC   5.248274  70.203751 25.648423 -6.454238  -5.834333  
GTEX-113JC  16.836411 -36.001744 10.849349  4.437152  -3.455231  
3)splicing: a matrix  
   ```splicing[1:5,1:10]```  
     #Chr start   end              ID GTEX-1LVAO GTEX-1AX9K GTEX-1GN73    GTEX-RM2N GTEX-111YS GTEX-1R9PN  
1 chr1 29552 29553 ENSG00000227232  1.5627810   1.617492  1.8133847  0.003239502  0.4766374  0.9620965  
2 chr1 29552 29553 ENSG00000227232 -1.3877237  -1.650869 -1.7964716 -0.004566534 -0.3512171 -0.8436864  
3 chr1 29552 29553 ENSG00000227232 -0.7163063  -2.124633  0.4795690 -0.005269084  1.7039249 -0.2537381  
4 chr1 29552 29553 ENSG00000227232  1.1117534   1.685548 -0.4241061 -0.014246458 -1.8304159  0.1309287  
5 chr1 29552 29553 ENSG00000227232  1.7357796   2.367227  0.6758256  0.156601906  0.2407400  1.1691889  

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
