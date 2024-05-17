simulate_blood_brain = function(){
    library(glmnet)
    load('exp_matrix.Rdata') #expression matrix
    blood_exp = exp.matrix[['Whole_Blood']]
    brain_exp = exp.matrix[['Brain_Caudate_basal_ganglia']]
    rm(exp.matrix)
    
    co_sams = intersect(row.names(blood_exp), row.names(brain_exp))
    blood_exp = blood_exp[co_sams, ]
    blood_exp_pc = prcomp(blood_exp)$x[, 1:60]
    brain_exp = as.data.frame(brain_exp[co_sams, ])
    brain_genes = colnames(brain_exp)
    rs = foreach(i = 1:length(brain_genes), .errorhandling = 'pass', .packages = 'glmnet') %dopar% {
        pred_gene = brain_genes[i]
        y = brain_exp[, pred_gene]
        model = cv.glmnet(
            x = as.matrix(blood_exp_pc), y = as.matrix(y), family = "gaussian", nfolds = 5, keep = T, alpha = 1
        )
        y_hat = model$fit.preval[, model$lambda == model$lambda.min]
        r = cor(y, y_hat)
        return(list('r' = r, 'beta' = model$glmnet.fit$beta[, model$lambda == model$lambda.min]))
    }
    
    rs_df = foreach(i = 1:length(brain_genes), .combine = 'rbind') %do% {
        return(data.frame('gene' = brain_genes[i], 'r' = rs[[i]][['r']]))
    }
    rs_df = rs_df[order(rs_df$r, decreasing = T), ]
    row.names(rs_df) = 1:dim(rs_df)[1]
    
    ### select 100 genes
    select_genes = NULL
    run = T
    gene_i = 1
    cor_threshold = 0.3
    while (run) {
        if (length(select_genes) >= 99 | rs_df[gene_i, 'r'] < cor_threshold){
            run = F
        }
        if (is.null(select_genes)){
            select_genes = rs_df$gene[gene_i]
        } else{
            current_gene = rs_df$gene[gene_i]
            pearson_r = c()
            for (before_gene in select_genes){
                pearson_r = c(pearson_r, cor(brain_exp[, before_gene], brain_exp[, current_gene]))
            }
            if (all(pearson_r < cor_threshold)){
                select_genes = c(select_genes, current_gene)
            }
        }
        gene_i = gene_i + 1
    }
    
    ### save beta
    select_gene_beta = NULL
    for (gene in select_genes){
        gene_beta = rs[[match(gene, brain_genes)]][['beta']]
        if (is.null(select_gene_beta)){
            select_gene_beta = gene_beta
        } else{
            select_gene_beta = rbind(select_gene_beta, gene_beta)
        }
    }
    row.names(select_gene_beta) = select_genes
    return(select_gene_beta)
}

simulate_brain_bmi = function(){
    select_gene_beta = simulate_blood_brain()
    select_genes = row.names(select_gene_beta)
    load('exp_matrix.Rdata') # expression matrix
    brain_exp = exp.matrix[['Brain_Caudate_basal_ganglia']]
    bmi = as.data.frame(data.table::fread('GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt'))[, c('SUBJID','BMI')] # BMI from GTEx
    
    co_sams = intersect(row.names(brain_exp), bmi$SUBJID)
    brain_exp = as.data.frame(brain_exp[co_sams, select_genes])
    bmi = bmi[match(co_sams, bmi$SUBJID), ]
    model = glmnet::cv.glmnet(
        x = as.matrix(brain_exp), y = as.matrix(bmi[, 'BMI']), family = "gaussian", nfolds = 5, keep = T, alpha = 1
    )
    y_hat = model$fit.preval[, model$lambda == model$lambda.min]
    print(paste0('brain -> bmi, r is ', cor(bmi[, 'BMI'], y_hat)))
    brain_bmi_beta = model$glmnet.fit$beta[, model$lambda == model$lambda.min]
    return(brain_bmi_beta)
}

simulate_dataset = function(num = 2){
    select_gene_beta = simulate_blood_brain()
    brain_bmi_beta = simulate_brain_bmi()
    simulation = list()
    for (i in  1:num){
        set.seed(413*i)
        rand_blood_exp = matrix(rnorm(300*60), nrow = 300, ncol = 60)
        simulate_brain = rand_blood_exp %*% t(select_gene_beta)
        colnames(simulate_brain) = row.names(select_gene_beta)
        simulate_bmi = simulate_brain %*%  matrix(brain_bmi_beta, ncol = 1)
        simulation[[i]] = list(
            'rand_blood_exp' = rand_blood_exp, 'simulate_brain' = simulate_brain, 'simulate_bmi' = simulate_bmi
        )
    }
    return(simulation)
}

simulate_validation = function(){
    rand_dataset = simulate_dataset()
    ### build prediction on dataset 1
    train_blood_exp = rand_dataset[[1]][['rand_blood_exp']]
    train_brain_exp = rand_dataset[[1]][['simulate_brain']]
    train_bmi = rand_dataset[[1]][['simulate_bmi']]
    pred_genes = colnames(train_brain_exp)
    train_pred_brain_exp = NULL
    blood_brain_pred_beta = NULL
    for (gene in pred_genes){
        model = glmnet::cv.glmnet(
            x = as.matrix(train_blood_exp), y = as.matrix(train_brain_exp[, gene]), family = "gaussian", nfolds = 5, keep = T, alpha = 1
        )
        y_hat = as.matrix(model$fit.preval[, model$lambda == model$lambda.min])
        model_beta = as.matrix(model$glmnet.fit$beta[, model$lambda == model$lambda.min])
        if (is.null(train_pred_brain_exp)){
            train_pred_brain_exp = y_hat
            blood_brain_pred_beta = model_beta
        } else{
            train_pred_brain_exp = cbind(train_pred_brain_exp, y_hat)
            blood_brain_pred_beta = cbind(blood_brain_pred_beta, model_beta)
        }
    }
    colnames(train_pred_brain_exp) = pred_genes
    
    ### DGE analysis on prediction of dataset 1
    dge_rs = data.frame()
    for (gene in pred_genes){
        dge_data = data.frame(
            'y' = train_bmi[, 1], 'x' = train_pred_brain_exp[, gene]
        )
        dge_lm = lm(y~x, data = dge_data)
        dge_rs[gene, 'train.beta'] = summary(dge_lm)$coefficients['x', 'Estimate']
    }
    
    ### validate on dataset 2
    val_blood_exp = rand_dataset[[2]][['rand_blood_exp']]
    val_brain_exp = rand_dataset[[2]][['simulate_brain']]
    val_bmi = rand_dataset[[2]][['simulate_bmi']]
    val_pred_brain_exp = val_blood_exp %*% blood_brain_pred_beta
    colnames(val_pred_brain_exp) = pred_genes
    for (gene in pred_genes){
        dge_data = data.frame(
            'y' = val_bmi[, 1], 'x' = val_pred_brain_exp[, gene]
        )
        dge_lm = lm(y~x, data = dge_data)
        dge_rs[gene, 'val.beta'] = summary(dge_lm)$coefficients['x', 'Estimate']
    }
    
    ### compare 2
    library(ggplot2)
    gg = ggplot(data = dge_rs, mapping = aes(x = train.beta, y = val.beta)) +
        geom_point(color = '#3e606f') +
        geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
        xlab('DGE-coefficients of pseudo GTEx') + ylab('DGE-coefficients of pseudo external dataset') +
        theme_bw() + coord_fixed() +
        theme(
            plot.margin = margin(20, 20, 20, 20)
        )
    gg
    pdf('/data/g_gamazon_lab/zhoud2/cross_tissue_prediction/revision/output/figures/simulation_val.pdf', width = 6, height = 6)
    plot(gg)
    dev.off()
    gg
}

