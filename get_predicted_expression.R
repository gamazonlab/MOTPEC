####get predicted expression using model trained by us
library(foreach)
library(doMC)
registerDoMC(cores=max(detectCores() - 1, 1))

##optional input data
splicing_index = T #whether add splicing
APA_index = T #whether add APA event
eSNP_index = T #whether add genetic variants

load('~/MOTPEC/data/output/all_model_beta.Rdata') #get beta of lasso
tissue.names = names(all_model_beta)
load('~/MOTPEC/data/iutput/WB.Rdata') #WB, partial infomation of WB
load('~/MOTPEC/data/iutput/sam_demo.Rdata') #dummy.demo.info, demography variables
load('~/MOTPEC/data/iutput/exp_matrix.Rdata') #exp.matrix, all tissues' expression
blood.exp <- exp.matrix[["Whole_Blood"]]
rm(exp.matrix)

args = as.numeric(commandArgs(TRUE))
tissue = tissue.names[args]
tissue.beta = all_model_beta[[tissue]]
pred.genes = names(tissue.beta)
blood.samples <- row.names(blood.exp)
demo.samples = row.names(dummy.demo.info)

samples <- intersect(blood.samples, demo.samples)
if (splicing_index){
  load('~/MOTPEC/data/iutput/splicing.Rdata') #splicing
  loci_path <- '~/MOTPEC/data/iutput/gencode.v32.GRCh37.txt' #loci profiles
  loci <- read.csv(loci_path, sep = '\t') #loci infomation
} 
if (APA_index){
  load('~/MOTPEC/data/iutput/apa_pc_data.Rdata') #apa 20 pcs
  apa.samples = row.names(apa_pc_data)
  samples = intersect(apa.samples, samples)
  apa.pc = apa_pc_data[samples, ]
} 
if (eSNP_index){
  loci_path <- '~/MOTPEC/data/iutput/gencode.v32.GRCh37.txt' #loci profiles
  loci <- read.csv(loci_path, sep = '\t') #loci infomation
}

co.blood.exp <- blood.exp[samples, ]
blood.pc <- WB[['blood.pca']][samples, ]
blood.tf.pc <- WB[['blood.tf.pca']][samples, ]

rs = foreach(i = 1:length(pred.genes), .combine = 'cbind') %dopar% {
  gene = pred.genes[i]
  model.beta = tissue.beta[[gene]]
  target.gene = co.blood.exp[, gene]
  module <- WB[['colors']][[gene]]
  module.pca <- WB[['modules_pca']][[module]][samples, ]
  x = cbind(target.gene, dummy.demo.info[samples, ], blood.pc, blood.tf.pc, module.pca)
  if (eSNP_index){
    local.snp <- getGeneSnpPCA(loci, gene, 10^6)
    local.snp <- local.snp[samples, ]
    x = cbind(x, local.snp)
  } else{
    local.snp = matrix(data = 0, nrow = length(samples), ncol = 20)
    x = cbind(x, local.snp)
  }
  if (splicing_index){
    splice.matrix <- splicing[ , c(colnames(splicing)[1:4], samples)]
    local.splice <- getGeneSplicePCA(loci, splice.matrix, gene, 10^5)
    x = cbind(x, local.splice)
  } else{
    local.splice = matrix(data = 0, nrow = length(samples), ncol = (nrow(model.beta) - 86))
    x = cbind(x, local.splice)
  }
  if (APA_index){
    x = cbind(x, apa.pc)
  } else{
    apa.pc = matrix(data = 0, nrow = length(samples), ncol = 20)
    x = cbind(x, apa.pc)
  }
  
  y = as.matrix(x) %*% as.matrix(model.beta)
  colnames(y) = gene
  return(y)
}

save(rs, file = paste0('~/MOTPEC/data/output/tissue_specific_exp/', tissue, '.Rdata'))
