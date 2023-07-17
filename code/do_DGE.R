work_dir = '~/MOTPEC'

load('~/MOTPEC/data/output/exp_pred.Rdata') #MOTPEC prediction
load('~/MOTPEC/data/iutput/exp_matrix.Rdata') #observed expression
load('~/MOTPEC/data/iutput/bmi.Rdata') #bmi phenotype
load('~/MOTPEC/data/iutput/sam_demo.Rdata') #demography variables

library(foreach)
library(doMC)
registerDoMC(cores=max(detectCores() - 1, 1))

tissue.names <- names(exp_pred)

all_tissues_beta <- list()
blood <- exp.matrix[['Whole_Blood']]
for (tissue in tissue.names){
  pred.exp <- exp_pred[[tissue]]
  pred.genes <- colnames(pred.exp) #the gene which predicted must in blood and observe
  pred.samples <- row.names(pred.exp)
  bmi.samples <- row.names(bmi)
  demo.samples <- row.names(dummy.demo.info)
  samples <- intersect(intersect(pred.samples, bmi.samples), demo.samples)
  
  co.blood <- blood[samples, ]
  co.tissue <- exp.matrix[[tissue]][samples, ]
  co.pred.tissue <- pred.exp[samples, ]
  co.bmi <- bmi[samples, 'BMI']
  co.demo <- dummy.demo.info[samples, ]
  
  rs <- foreach(i = 1:length(pred.genes), .combine = 'rbind') %dopar% {
    rs_ <- data.frame()
    gene <- pred.genes[i]
    
    observe.data <- cbind(co.bmi, co.tissue[, gene], co.demo)
    colnames(observe.data) <- c('y', 'x1', 'x2', 'x3', 'x41', 'x42', 'x43')
    observe.lm <- lm(y~., data = observe.data)
    observe.sum <- summary(observe.lm)
    rs_[gene, 'observe.beta'] <- observe.sum$coefficients['x1', 'Estimate']
    rs_[gene, 'observe.se'] <- observe.sum$coefficients['x1', 'Std. Error']
    rs_[gene, 'observe.p'] <- observe.sum$coefficients['x1', 'Pr(>|t|)']
    
    pred.data <- cbind(co.bmi, co.pred.tissue[, gene], co.demo)
    colnames(pred.data) <- c('y', 'x1', 'x2', 'x3', 'x41', 'x42', 'x43')
    pred.lm <- lm(y~., data = pred.data)
    pred.sum <- summary(pred.lm)
    rs_[gene, 'pred.beta'] <- pred.sum$coefficients['x1', 'Estimate']
    rs_[gene, 'pred.se'] <- pred.sum$coefficients['x1', 'Std. Error']
    rs_[gene, 'pred.p'] <- pred.sum$coefficients['x1', 'Pr(>|t|)']
    
    blood.data <- cbind(co.bmi, co.blood[, gene], co.demo)
    colnames(blood.data) <- c('y', 'x1', 'x2', 'x3', 'x41', 'x42', 'x43')
    blood.lm <- lm(y~., data = blood.data)
    blood.sum <- summary(blood.lm)
    rs_[gene, 'blood.beta'] <- blood.sum$coefficients['x1', 'Estimate']
    rs_[gene, 'blood.se'] <- blood.sum$coefficients['x1', 'Std. Error']
    rs_[gene, 'blood.p'] <- blood.sum$coefficients['x1', 'Pr(>|t|)']
    return(rs_)
  }
  all_tissues_beta[[tissue]] <- rs
}
save(all_tissues_beta, file = '~/MOTPEC/data/output/all_tissues_beta.Rdata')

#get whole blood beta = b.beta
bmi.samples <- row.names(bmi)
demo.samples <- row.names(dummy.demo.info)
samples <- intersect(intersect(row.names(blood), bmi.samples), demo.samples)
blood_genes = colnames(blood)
rs_ = data.frame()
for (i in 1:length(blood_genes)){
  gene = blood_genes[i]
  co.bmi <- bmi[samples, 'BMI']
  co.demo <- dummy.demo.info[samples, ]
  x = cbind(co.bmi, blood[samples, gene], co.demo)
  colnames(x) = c('y', 'x1', 'x2', 'x3', 'x41', 'x42', 'x43')
  blood.lm <- lm(y~., data = x)
  blood.sum <- summary(blood.lm)
  rs_[gene, 'blood.beta'] <- blood.sum$coefficients['x1', 'Estimate']
  rs_[gene, 'blood.se'] <- blood.sum$coefficients['x1', 'Std. Error']
  rs_[gene, 'blood.p'] <- blood.sum$coefficients['x1', 'Pr(>|t|)']
}
b.beta = rs_
save(b.beta, file = '~/MOTPEC/data/output/b_beta.Rdata')
