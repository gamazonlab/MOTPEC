####get predicted expression using model trained by us
library(foreach)
library(doMC)
library(optparse)
registerDoMC(cores=max(detectCores() - 1, 1))

option_list = list(
  make_option('--input_dir', action="store", default=NA, type='character'),
  make_option('--output_dir', action="store", default=NA, type='character'),
  make_option('--demo_index', action = 'store', default = F, type = 'logical'),
  make_option('--APA_index', action = 'store', default = F, type = 'logical'),
  make_option('--eSNP_index', action = 'store', default = F, type = 'logical'),
  make_option('--plink', action="store", default=NA, type='character'),
  make_option('--plink_wd', action="store", default=NA, type='character'),
  make_option('--plink_all_sam_name', action="store", default=NA, type='character'),
  make_option('--splicing_index', action = 'store', default = F, type = 'logical'),
  make_option('--args', action = 'store', default = NA, type = 'integer')
)

opt = parse_args(OptionParser(option_list=option_list))
print('Your arguments are as followes:')
print(opt)

##get arguments
input_dir = opt$input_dir
output_dir = opt$output_dir
demo_index = opt$demo_index
APA_index = opt$APA_index
eSNP_index = opt$eSNP_index
plink = opt$plink
plink_wd = opt$plink_wd
plink_all_sam_name = opt$plink_all_sam_name
splicing_index = opt$splicing_index
args = opt$args

if(!dir.exists(input_dir)){
  dir.create(input_dir)
}
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

load(paste0(input_dir, '/all_model_beta.Rdata')) #get beta of lasso, we provide
load(paste0(input_dir, '/co_exp_net_and_blood_pc.Rdata')) #get rotation of pc, we provide
tissue.names = names(all_model_beta)
load(paste0(input_dir, '/blood_exp.Rdata')) #blood.exp

tissue = tissue.names[args]
tissue.beta = all_model_beta[[tissue]]
pred.genes = names(tissue.beta)
blood.samples <- row.names(blood.exp)

samples <- blood.samples
if (demo_index){
  load(paste0(input_dir, '/sam_demo.Rdata')) #dummy.demo.info, demography variables
  demo.samples = row.names(dummy.demo.info)
  samples = intersect(samples, demo.samples)
}
if (splicing_index){
  #a function for getting scale and chromosome of gene
  getGeneLeftRight <- function(loci, gene_id, distance){
    row_info <- loci[which(loci$geneid == gene_id),]
    chr <- as.character(row_info[1, 'chr'])
    left <- max(row_info[1, 'left'] - distance, min(loci[which(loci$chr == chr),]$left))
    right <- min(row_info[1, 'right'] + distance, max(loci[which(loci$chr == chr),]$right))
    return(c(chr, left, right))
  }
  
  #a function for getting PCs of local splicing
  getGeneSplicePCA <- function(loci, splicing, gene_id, distance){
    loc.info <- getGeneLeftRight(loci, gene_id, distance)
    gene.chr <- as.character(loc.info[1])
    gene.left <- as.numeric(loc.info[2])
    gene.right <- as.numeric(loc.info[3])
    temp.loci <- loci[which(loci[,1] == gene.chr),]
    temp.loci <- temp.loci[which((temp.loci[,2] > gene.left) & (temp.loci[,3] < gene.right)),]
    gene.scale <- temp.loci[,6] #在范围内的基因id
    
    splice.gene.ids <- names(table(splicing[,4]))
    gene.comm <- intersect(gene.scale, splice.gene.ids)
    splice.dt <- NULL #splice df
    for (gene_ in gene.comm){
      temp.splice.dt <- splicing[which(splicing[,4] == gene_),]
      splice.dt <- rbind(splice.dt, temp.splice.dt)
    }
    splice.dt <- t(splice.dt[,-c(1:4)])
    
    if (length(colnames(splice.dt)) > 20){
      splice.dt <- data.frame(prcomp(splice.dt, rank = 20)$x)
    }
    
    return(splice.dt)
  }
  
  load(paste0(input_dir, '/splicing.Rdata')) #splicing
  loci_path <- paste0(input_dir, '/gencode.v32.GRCh37.txt') #loci profiles
  loci <- read.csv(loci_path, sep = '\t') #loci infomation
} 

if (APA_index){
  load(paste0(input_dir, '/apa_pc_data.Rdata')) #apa 20 pcs
  apa.samples = row.names(apa_pc_data)
  samples = intersect(apa.samples, samples)
  apa.pc = apa_pc_data[samples, ]
} 

if (eSNP_index){
  #a function for getting scale and chromosome of gene
  getGeneLeftRight <- function(loci, gene_id, distance){
    row_info <- loci[which(loci$geneid == gene_id),]
    chr <- as.character(row_info[1, 'chr'])
    left <- max(row_info[1, 'left'] - distance, min(loci[which(loci$chr == chr),]$left))
    right <- min(row_info[1, 'right'] + distance, max(loci[which(loci$chr == chr),]$right))
    return(c(chr, left, right))
  }
  
  #a function for getting PCs of local eSNP
  getGeneSnpPCA <- function(loci, gene_id, distance){
    rs <- getGeneLeftRight(loci, gene_id, distance)
    gene_id <- substr(gene_id, 1, 15)
    if (is.na(rs[1])){
      return(NULL)
    } else{
      #get basic info
      chr <- substr(rs[1], 4, nchar(rs[1]))
      left <- rs[2]
      right <- rs[3]
      #plink <- '/gpfs52/home/zhoud2/tools/plink/plink'
      #plink_wd = '/gpfs52/data/g_gamazon_lab/zhoud2/cross_tissue_prediction/genotype/'
      #plink_all_sam_name = 'gtex_v8_all'
      plink_str_snp <- paste0(plink, ' --bfile ', plink_wd, plink_all_sam_name,' --chr ', chr, 
                              ' --from-bp ', left, ' --to-bp ',
                              right, ' --make-bed --out ', plink_wd, gene_id)
      plink_str_snp_pca <- paste0(plink, 
                                  ' --bfile ', plink_wd, gene_id, ' --pca 20 --out ', plink_wd,
                                  paste0(gene_id, '_pca'))
      system(plink_str_snp, wait = T)
      system(plink_str_snp_pca, wait = T)
      #get pca
      snp_pca <- read.csv(paste0(plink_wd, paste0(gene_id, '_pca'), '.eigenvec'),
                          sep = ' ', header = F, row.names = 1)[,-1]
      #drop file
      drop.file.name <- c('.bed', '.bim', '.fam', '.log', '.nosex')
      drop.file.pca.name <- c('.eigenval', '.eigenvec', '.log', '.nosex')
      for (i in drop.file.name) {
        drop.str <- paste0('rm -f ', plink_wd, gene_id, i)
        system(drop.str, wait = T)
      }
      for (i in drop.file.pca.name){
        drop.str <- paste0('rm -f ', plink_wd, gene_id, '_pca', i)
        system(drop.str, wait = T)
      }
      return(snp_pca)
    }
  }
  loci_path <- paste0(input_dir, '/gencode.v32.GRCh37.txt') #loci profiles
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
  tryCatch({
    print(paste0('Now we are predicting the No.', i, ' gene of ', tissue, '!'))
    gene = pred.genes[i]
    model.beta = tissue.beta[[gene]]
    target.gene = as.matrix(co.blood.exp[, gene])
    
    #use module profiles we provide to calculate the module PCA
    module <- co_exp_net_and_blood_pc[['colors']][[gene]] #module of gene
    module.pca.rotation <- co_exp_net_and_blood_pc[['modules_pca']][[module]]
    module.genes = names(co_exp_net_and_blood_pc[['modules']][[module]])
    pc.co.module.genes = intersect(colnames(co.blood.exp), module.genes)
    module.pca = (as.matrix(co.blood.exp[,pc.co.module.genes]) %*% as.matrix(module.pca.rotation[pc.co.module.genes,]))[,1:20]
    
    #NOTE: the order of input features is target.gene,dummy.demo.info[samples, ],apa.pc,blood.pc,blood.tf.pc,module.pca,local.snp,local.splice
    #because the beta we provide is in this order, 
    #if you use your beta, please make sure that features' order is same with beta
    x = target.gene
    if (demo_index){
      x = cbind(x, as.matrix(dummy.demo.info[samples, ]))
    } else{
      x = cbind(x, matrix(data = 0, nrow = length(samples), ncol = 5))
    }
    
    if (APA_index){
      x = cbind(x, as.matrix(apa.pc))
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
    #if (dim(y)[1] != 200){y = data.frame(matrix(NA, nrow = 200))}
    colnames(y) = gene
    #row.names(y) = 1:length(samples)
    return(y)
  }, error = function(e){
    #print(e)
  })
  
}

write.table(rs, file = paste0(output_dir, '/', tissue, '.txt'), sep = '\t')
print('Predicting is OK!')
