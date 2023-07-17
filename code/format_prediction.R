###function : get rs from prediction results
get_model_rs = function(rs){
  tissue_rs = rs
  pcc_df = data.frame()
  pred_exp = data.frame()
  model_beta = list()
  for (i in 1:(length(tissue_rs) - 2)){
    if (is.null(tissue_rs[[i]])){
      next
    } else{
      #print(i)
      gene_name = names(tissue_rs[[i]])
      gene_rs = tissue_rs[[i]][[1]]
      pcc_df[gene_name, 'baseline'] = gene_rs[['baseline.r']]
      pcc_df[gene_name, 'MOTPEC'] = gene_rs[['MOTPEC.r']]
      model_beta[[gene_name]] = gene_rs[['MOTPEC.beta']]
      if (dim(pred_exp)[1] == 0){
        colnames(gene_rs[['MOTPEC.prediction']]) = gene_name
        pred_exp = gene_rs[['MOTPEC.prediction']]
      } else{
        colnames(gene_rs[['MOTPEC.prediction']]) = gene_name
        pred_exp = cbind(pred_exp, gene_rs[['MOTPEC.prediction']])
      }
    }
  }
  sample_num = tissue_rs[[length(tissue_rs)]]
  model_time = tissue_rs[[length(tissue_rs) - 1]]
  
  info = list()
  info[['pcc_df']] = pcc_df
  info[['pred_exp']] = pred_exp
  info[['model_beta']] = model_beta
  info[['sample_num']] = sample_num
  info[['pred_gene_num']] = dim(pred_exp)[2]
  info[['gene_num']] = length(tissue_rs) - 2
  info[['model_time']] = model_time
  return(info)
}

###get pcc, pred_exp, beta, sam, gene, time of each tissue
file_names = dir('~/MOTPEC/data/output/prediction_model_rs/')
pcc_2models = list()
exp_pred = list()
all_model_beta = list()
sam_gene_time = data.frame()
for (file_name in file_names){
  tissue = substr(file_name, 1, nchar(file_name) - 6)
  print(tissue)
  load(paste0('~/MOTPEC/data/output/prediction_model_rs/', file_name))
  temp_rs = get_model_rs(rs)
  rm(rs)
  pcc_2models[[tissue]] = temp_rs[['pcc_df']]
  exp_pred[[tissue]] = temp_rs[['pred_exp']]
  all_model_beta[[tissue]] = temp_rs[['model_beta']]
  sam_gene_time[tissue, 'sample_num'] = temp_rs[['sample_num']]
  sam_gene_time[tissue, 'pred_gene_num'] = temp_rs[['pred_gene_num']]
  sam_gene_time[tissue, 'gene_num'] = temp_rs[['gene_num']]
  sam_gene_time[tissue, 'model_time'] = temp_rs[['model_time']]
  rm(temp_rs)
}
sam_gene_time$pred_ratio = sam_gene_time$pred_gene_num / sam_gene_time$gene_num

###save these results
save(pcc_2models, file = '~/MOTPEC/data/output/pcc_2models.Rdata')
save(exp_pred, file = '~/MOTPEC/data/output/exp_pred.Rdata')
save(all_model_beta, file = '~/MOTPEC/data/output/all_model_beta.Rdata')
save(sam_gene_time, file = '~/MOTPEC/data/output/sam_gene_time.Rdata')
