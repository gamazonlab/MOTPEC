#a function for getting scale and chromosome of gene
#loci : a dataframe for gene location profiles
#gene_id : geneid for target gene
#return : vector, c(chr, left, right)
getGeneLeftRight <- function(loci, gene_id, distance){
  row_info <- loci[which(loci$geneid == gene_id),]
  chr <- as.character(row_info[1, 'chr'])
  left <- max(row_info[1, 'left'] - distance, min(loci[which(loci$chr == chr),]$left))
  right <- min(row_info[1, 'right'] + distance, max(loci[which(loci$chr == chr),]$right))
  return(c(chr, left, right))
}

#a function for getting PCs of local eSNP
#loci : a dataframe for gene location profiles
#gene_id : geneid for target gene
#return : dataframe, PC(local_snp)
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
    plink <- '/gpfs52/home/zhoud2/tools/plink/plink'
    plink_str_snp <- paste0(plink, ' --bfile genotype/gtex_v8_all --chr ', chr, 
                            ' --from-bp ', left, ' --to-bp ',
                            right, ' --make-bed --out genotype/', gene_id)
    plink_str_snp_pca <- paste0(plink, 
                                ' --bfile genotype/', gene_id, ' --pca 20 --out genotype/',
                                paste0(gene_id, '_pca'))
    system(plink_str_snp, wait = T)
    system(plink_str_snp_pca, wait = T)
    #get pca
    snp_pca <- read.csv(paste0('genotype/', paste0(gene_id, '_pca'), '.eigenvec'),
                        sep = ' ', header = F, row.names = 1)[,-1]
    #drop file
    drop.file.name <- c('.bed', '.bim', '.fam', '.log', '.nosex')
    drop.file.pca.name <- c('.eigenval', '.eigenvec', '.log', '.nosex')
    for (i in drop.file.name) {
      drop.str <- paste0('rm -f genotype/', gene_id, i)
      system(drop.str, wait = T)
    }
    for (i in drop.file.pca.name){
      drop.str <- paste0('rm -f genotype/', gene_id, '_pca', i)
      system(drop.str, wait = T)
    }
    return(snp_pca)
  }
}

#a function for getting PCs of local splicing
#loci : a dataframe for gene location profiles
#splicing : a dataframe of splicing profiles
#gene_id : geneid for target gene
#distance : distance for local profiles' scale
#return : splice pca
getGeneSplicePCA <- function(loci, splicing, gene_id, distance){
  loc.info <- getGeneLeftRight(loci, gene_id, distance)
  gene.chr <- as.character(loc.info[1])
  gene.left <- as.numeric(loc.info[2])
  gene.right <- as.numeric(loc.info[3])
  temp.loci <- loci[which(loci[,1] == gene.chr),]
  temp.loci <- temp.loci[which((temp.loci[,2] > gene.left) & (temp.loci[,3] < gene.right)),]
  gene.scale <- temp.loci[,6] #gene id in scale
  
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
