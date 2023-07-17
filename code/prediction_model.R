work_dir = '~/MOTPEC'

source("~/MOTPEC/code/functions.R")
library(glmnet)
library(data.table)
library(caret)
library(foreach)
library(doMC)
registerDoMC(cores=max(detectCores() - 1, 1))

##optional input data
splicing_index = T #whether add splicing
APA_index = T #whether add APA event
eSNP_index = T #whether add genetic variants

load('~/MOTPEC/data/iutput/exp_matrix.Rdata') #exp.matrix, all tissues' expression
load('~/MOTPEC/data/iutput/WB.Rdata') #WB, partial infomation of WB
load('~/MOTPEC/data/iutput/sam_demo.Rdata') #dummy.demo.info, demography variables

tissue.names <- names(exp.matrix)
blood.exp <- exp.matrix[["Whole_Blood"]]
tissue.names <- tissue.names[-49] #del whole blood
args = as.numeric(commandArgs(TRUE))
tissue = tissue.names[args]
tissue.exp <- exp.matrix[[tissue]]
tissue.samples <- row.names(tissue.exp)
blood.samples <- row.names(blood.exp)
demo.samples = row.names(dummy.demo.info)

samples <- intersect(tissue.samples, intersect(blood.samples, demo.samples))
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

co.tissue.exp <- tissue.exp[samples, ]
co.blood.exp <- blood.exp[samples, ]
blood.pc <- WB[['blood.pca']][samples, ]
blood.tf.pc <- WB[['blood.tf.pca']][samples, ]
pred.genes <- colnames(co.tissue.exp)

start = Sys.time()

rs <- foreach(i = 1:length(pred.genes)) %dopar% {
  tryCatch({
    gene <- pred.genes[i]
    results = data.frame()
    seed_num <- as.integer(substr(gene, 5, nchar(gene)))
    
    #baseline model(null model)
    target.gene <- co.blood.exp[, gene]
    y <- co.tissue.exp[, gene]
    baseline.r <- cor(target.gene, y)
    results[[gene]][['baseline.r']] = baseline.r
    
    #MOTPEC prediction model
    module <- WB[['colors']][[gene]]
    module.pca <- WB[['modules_pca']][[module]][samples, ]
    x = cbind(target.gene, dummy.demo.info[samples, ], blood.pc, blood.tf.pc, module.pca)
    if (eSNP_index){
      local.snp <- getGeneSnpPCA(loci, gene, 10^6)
      local.snp <- local.snp[samples, ]
      x = cbind(x, local.snp)
    }
    if (splicing_index){
      splice.matrix <- splicing[ , c(colnames(splicing)[1:4], samples)]
      local.splice <- getGeneSplicePCA(loci, splice.matrix, gene, 10^5)
      x = cbind(x, local.splice)
    }
    if (APA_index){
      x = cbind(x, apa.pc)
    }
    
    #set penalty factors
    pfs <- rep(1, length(colnames(x)))
    pfs[1] <- 0
    
    set.seed(seed_num)
    ans_our_model = cv.glmnet(x = as.matrix(x), y = as.matrix(y),
                              penalty.factor = pfs,
                              family = "gaussian", 
                              nfolds = 5,
                              keep = T,
                              alpha = 1)
    
    y_hat = ans_our_model$fit.preval[,which(ans_our_model$lambda == ans_our_model$lambda.min)] #predicted value
    performance <- cor.test(y, y_hat) #pearson's r
    y_hat = as.data.frame(y_hat)
    row.names(y_hat) = samples
    
    results[[gene]][['MOTPEC.r']] = as.numeric(performance$estimate[[1]])
    results[[gene]][['MOTPEC.prediction']] = y_hat
    results[[gene]][['MOTPEC.beta']] = as.data.frame(as.numeric(ans_our_model$glmnet.fit$beta[, which(ans_our_model$lambda == ans_our_model$lambda.min)]))
    
    return(results)
  }, error = function(e){
    #print(e)
  })
}

end = Sys.time()
rs[['time']] = end - start
rs[['sample_num']] = length(samples)
save(rs, file = paste0('~/MOTPEC/data/output/prediction_model_rs/', tissue, '.Rdata'))
