ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

cd /gpfs52/home/zhoud2/script/xuyue/cross_tissue_pred_exp/MOTPEC/code/get_predicted_expression/

Rscript get_predicted_expression.R \
  --input_dir '/gpfs52/data/g_gamazon_lab/zhoud2/cross_tissue_prediction/test_MOTPEC/' \
  --output_dir '/gpfs52/data/g_gamazon_lab/zhoud2/cross_tissue_prediction/test_MOTPEC/output/' \
  --args 1 \
  --demo_index T \
  #--APA_index T \
  #--eSNP_index T \
  #--plink '/gpfs52/home/zhoud2/tools/plink/plink' \
  #--plink_wd '/gpfs52/data/g_gamazon_lab/zhoud2/cross_tissue_prediction/genotype/' \
  #--plink_all_sam_name 'gtex_v8_all' \
  #--splicing_index T 
  
