####get predicted expression using model trained by us
library(foreach)
library(doMC)
registerDoMC(cores=max(detectCores() - 1, 1))

##optional input data
splicing_index = T #whether add splicing
APA_index = T #whether add APA event
eSNP_index = T #whether add genetic variants

load('~/MOTPEC/data/provide/all_model_beta.Rdata') #get beta of lasso, we provide
load('~/MOTPEC/data/provide/co_exp_net_and_blood_pc.Rdata') #get rotation of pc, we provide
tissue.names = names(all_model_beta)
#dummy.demo.info and blood expression, you must prepare by yourself
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

###use loading score we provide to calculate the PCA
co.blood.exp <- blood.exp[samples, ]
blood.pca.rotation = co_exp_net_and_blood_pc[['blood.pca.rotation']]
pc.co.genes = intersect(colnames(co.blood.exp), row.names(blood.pca.rotation))
blood.pc <- (as.matrix(co.blood.exp[,pc.co.genes]) %*% as.matrix(blood.pca.rotation[pc.co.genes,]))[,1:10]

blood.tf.pca.rotation = co_exp_net_and_blood_pc[['blood.tf.pca.rotation']]
pc.co.tf.genes = intersect(colnames(co.blood.exp), row.names(blood.tf.pca.rotation))
blood.tf.pc <- (as.matrix(co.blood.exp[,pc.co.tf.genes]) %*% as.matrix(blood.tf.pca.rotation[pc.co.tf.genes,]))[,1:10]

rs = foreach(i = 1:length(pred.genes), .combine = 'cbind') %dopar% {
  gene = pred.genes[i]
  model.beta = tissue.beta[[gene]]
  target.gene = co.blood.exp[, gene]
  
  #use module profiles we provide to calculate the module PCA
  module <- co_exp_net_and_blood_pc[['colors']][[gene]] #module of gene
  module.pca.rotation <- co_exp_net_and_blood_pc[['modules_pca']][[module]]
  module.genes = names(co_exp_net_and_blood_pc[['modules']])
  pc.co.module.genes = intersect(colnames(co.blood.exp), module.genes)
  module.pca = (as.matrix(co.blood.exp[,pc.co.module.genes]) %*% as.matrix(module.pca.rotation[pc.co.module.genes,]))[,1:20]
  
  #NOTE: the order of input features is target.gene,dummy.demo.info[samples, ],apa.pc,blood.pc,blood.tf.pc,module.pca,local.snp,local.splice
  #because the beta we provide is in this order, 
  #if you use your beta, please make sure that features' order is same with beta
  x = cbind(target.gene, dummy.demo.info[samples, ])
  if (APA_index){
    x = cbind(x, apa.pc)
  } else{
    apa.pc = matrix(data = 0, nrow = length(samples), ncol = 20)
    x = cbind(x, apa.pc)
  }
  x = cbind(x, blood.pc, blood.tf.pc, module.pca)
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
  
  y = as.matrix(x) %*% as.matrix(model.beta)
  colnames(y) = gene
  return(y)
}

save(rs, file = paste0('~/MOTPEC/data/output/tissue_specific_exp/', tissue, '.Rdata'))
