#get tissue names
args = as.numeric(commandArgs(TRUE))
load('~/MOTPEC/data/iutput/exp_matrix.Rdata')
tissue.names <- names(exp.matrix)
tissue.names <- tissue.names[args]
remove(exp.matrix)

src_dir <- '~/MOTPEC/code/do_TWAS/'
db_path <- '~/MOTPEC/data/raw_data/GTEx_v8/db/'
cov_path <- '~/MOTPEC/data/raw_data/cov/'
gwas_path <- '~/MOTPEC/data/raw_data/giant'
gwas_name <- 'newbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz'
asso_out_path <- '~/MOTPEC/data/output/twas_rs'
for (tissue in tissue.names){
  db_name <- paste0('combo_', tissue, '.db')
  cov_name <- paste0('combo_', tissue, '.txt.gz')
  asso_out_name <- paste0(tissue, '.txt')
  
  run_str <- paste0('Rscript ', src_dir, '/predixcan_r.r ', 
                    '--asso_test ',
                    '--db_path ', db_path, '/', db_name, ' ',
                    '--cov_path ', cov_path, '/', cov_name, ' ',
                    '--gwas_path ', gwas_path, '/', gwas_name, ' ',
                    '--gwas_variant_col SNP --gwas_beta_col BETA --gwas_se_col SE --gwas_eff_allele_col Tested_Allele --gwas_ref_allele_col Other_Allele ',
                    '--asso_out_path ', asso_out_path, '/', asso_out_name, ' ',
                    '--parallel'
  )
  system(run_str, wait = T)
}
