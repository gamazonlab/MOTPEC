work_dir = '~/MOTPEC'

library(data.table)
library(WGCNA)
library(ggplot2)
library(readr)

### get normal exp
exp_dir <- '~/MOTPEC/data/raw_data/normal_exp'
exp.filename <- dir(exp_dir) #filename
exp.matrix <- list()

for (name in exp.filename){
  dir_ <- paste0(exp_dir, name)
  data_ <- as.data.frame(fread(dir_))
  for (i in 1:length(data_$gene_id)){
    data_$gene_id[i] <- strsplit(data_$gene_id[i], '[.]')[[1]][1]
  }
  row.names(data_) <- data_$gene_id
  data_ <- t(data_[ , -c(1:4)])
  tissue.name <- strsplit(name, '[.]')[[1]][1]
  exp.matrix[[tissue.name]] <- data_
}

exp.matrix.str <- '~/MOTPEC/data/iutput/exp_matrix.Rdata'
save(exp.matrix, file = exp.matrix.str)

### get partial whole blood profiles, including blood.pca, blood.tf.pca, modules_pca
blood <- exp.matrix[["Whole_Blood"]]
tf_path <- "MOTPEC/data/raw_data/TF.txt"
tf <- as.data.frame(fread(tf_path))
blood.genes <- colnames(blood)
blood.tfs <- intersect(blood.genes, tf[,2]) #全血转录因子
blood.pca <- prcomp(blood)$x[ , c(1:10)]
blood.tf.pca <- prcomp(blood[ , blood.tfs])$x[ , c(1:10)]

 ## do wgcna and get co-expression modules
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(blood, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.80,col="blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
abline(h=0.80, col="blue")
par(mfrow = c(1,1))
net = blockwiseModules(blood, power = 10, corType = "pearson", networkType = "unsigned", TOMType = "unsigned", minModuleSize = 30, numericLabels = T, mergeCutHeight = 0.25, verbose = 3)
colors <- net$colors
module_colors <- names(table(net$colors))
modules <- list()
for (module_color in module_colors){
  modules[[module_color]] <- which(colors == module_color)
}
modules_pca <- list()
for (module in module_colors){
  genes <- names(modules[[module]])
  x <- blood[,genes]
  modules_pca[[module]] <- prcomp(x)$x[,c(1:20)]
}

WB <- list()
WB[['blood.pca']] <- blood.pca
WB[['blood.tf.pca']] <- blood.tf.pca
WB[['colors']] <- colors
WB[['modules_pca']] <- modules_pca
WB[['modules']] <- modules
file.str <- '~/MOTPEC/data/iutput/WB.Rdata'
save(WB, file = file.str)

### get splicing
splicing_path <- 'MOTPEC/data/raw_data/Whole_Blood.v8.leafcutter_phenotypes.bed.gz'
splicing <- as.data.frame(fread(paste0(data_path, splicing_path))) #splicing infomation
for (i in 1:nrow(splicing)){
  splicing[i, 4] <- strsplit(strsplit(splicing[i, 4], ':')[[1]][5], '[.]')[[1]][1]
}
save(splicing, file = '~/MOTPEC/data/iutput/splicing.Rdata')

### get apa event
sam_attr = as.data.frame(fread('MOTPEC/data/raw_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'))
blood_sam_attr = sam_attr[sam_attr$SMTSD == 'Whole Blood', ]
raw_apa_data = as.data.frame(fread(unzip('MOTPEC/data/raw_data/Whole_Blood_All_PDUIs.zip', "Whole_Blood_All_PDUIs.txt")))
apa_data = t(raw_apa_data[,-1])
 ##count na
count_na = data.frame()
for (col in 1:ncol(apa_data)){
  nas = sum(is.na(apa_data[,col]))
  count_na[col, 'NA_num'] = nas
  count_na[col, 'NA_ratio'] = nas / (length(apa_data[,col]) - 1)
}
apa_data = apa_data[, which(count_na$NA_ratio < 0.2)]
for (i in 1:ncol(apa_data)){
  apa_data[, i][is.na(apa_data[, i])] = median(as.numeric(apa_data[, i]), na.rm = T)
}
apa_pc_data = as.data.frame(prcomp(apa_data)$x[, 1:20])
save(apa_pc_data, file = '~/MOTPEC/data/iutput/apa_pc_data.Rdata')

###get demography variables
demographic_path <- "~/MOTPEC/data/raw_data/sample_ga.txt"
demo.info <- as.data.frame(fread(demographic_path))[, c(1, 3, 4, 5)]
demo.info <- demo.info[-which(demo.info$RACE == 4 | demo.info$RACE == 99),]
demo.info$SEX = demo.info$SEX - 1
demo.info$RACE <- factor(demo.info$RACE)
dummy.demo <- dummyVars(' ~ RACE', data = demo.info)
dummy.demo.info <- as.data.frame(predict(dummy.demo, newdata = demo.info))
dummy.demo.info <- cbind(demo.info[, c(2, 3)], dummy.demo.info)
row.names(dummy.demo.info) <- demo.info[, 1]
save(dummy.demo.info, file = '~/MOTPEC/data/iutput/sam_demo.Rdata')

###get BMI phenotype
bmi_path <- '~/MOTPEC/data/raw_data/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt'
bmi <- as.data.frame(fread(bmi_path))[c('SUBJID','BMI')]
row.names(bmi) <- bmi[, 1]
save(bmi, file = '~/MOTPEC/data/iutput/bmi.Rdata')
