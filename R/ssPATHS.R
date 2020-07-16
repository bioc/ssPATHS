
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


# pathway activation score functions
library("dml")
library("data.table")
library("ROCR")
library("MESS")

get_hypoxia_genes <- function(){

    HIF_GENES = c("ENSG00000148926", "ENSG00000109107", "ENSG00000176171",
                  "ENSG00000104765", "ENSG00000074410", "ENSG00000107159",
                  "ENSG00000130635", "ENSG00000047457", "ENSG00000168209",
                  "ENSG00000129521", "ENSG00000111674", "ENSG00000104812",
                  "ENSG00000159399", "ENSG00000100292", "ENSG00000134333",
                  "ENSG00000113083", "ENSG00000123384", "ENSG00000185499",
                  "ENSG00000119950", "ENSG00000104419", "ENSG00000185633",
                  "ENSG00000124785", "ENSG00000152256", "ENSG00000114268",
                  "ENSG00000204531", "ENSG00000119938", "ENSG00000139832",
                  "ENSG00000141526", "ENSG00000117394", "ENSG00000103257",
                  "ENSG00000113739", "ENSG00000265972", "ENSG00000112715",
                  "ENSG00000186918", "ENSG00000117289")
    return(HIF_GENES)
}

get_gene_weights <- function(expression_matr, gene_ids){

    if(sum(colnames(expression_matr) %in% c("Y", "sample_id")) != 2){
        stop("Need column names Y and sample_id")
    }
    if(length(unique(expression_matr$Y)) != 2){
        stop("Y must be binary")
    }
    if(!(1 %in% expression_matr$Y) | !(0 %in% expression_matr$Y)){
        stop("Y must have both a 1 and 0 entry")
    }

    total_gene_ids = setdiff(colnames(expression_matr), c("Y", "sample_id"))

    # normalize
    dca_matr = expression_matr[,total_gene_ids]
    dca_data = as.matrix(t(scale(log10(t(dca_matr+1)))))
    dca_data[is.nan(dca_data)] <- 0
    row.names(dca_data) = expression_matr$sample_id
    genes_interest_idx = which(total_gene_ids %in% gene_ids)
    gene_ids_final = colnames(dca_data)[genes_interest_idx]

    dca_data = dca_data[,gene_ids_final]

    # set up DCA
    chunks = rep(1, nrow(expression_matr))
    chunks[expression_matr$Y!=0] = 2
    neglinks = matrix(c(0, 1, 1, 0), 2, 2)

    # get weights
    dca_res = dml::dca(data=dca_data, chunks=chunks, neglinks=neglinks)
    proj_vector = t(as.matrix(dca_res$DCA))

    # get projection
    dca_proj = dca_data %*% proj_vector
    dca_proj = data.frame(pathway_score=dca_proj, sample_id=row.names(dca_proj))
    dca_proj = merge(expression_matr[,c("sample_id", "Y")], dca_proj)
    dca_proj = dca_proj[order(dca_proj$pathway_score),]

    # check direction
    lower_score = mean(dca_proj$pathway_score[dca_proj$Y==0])
    upper_score = mean(dca_proj$pathway_score[dca_proj$Y==1])
    flipped = FALSE
    if(lower_score > upper_score){
        proj_vector = proj_vector * -1
        dca_proj = dca_data %*% proj_vector
        dca_proj = data.frame(pathway_score=dca_proj, sample_id=row.names(dca_proj))
        dca_proj = merge(expression_matr[,c("sample_id", "Y")], dca_proj)
        dca_proj = dca_proj[order(dca_proj$pathway_score),]
        flipped = TRUE
    }

    proj_vector_df = data.frame(gene_weight=proj_vector, gene_ids=gene_ids_final)

    return(list(proj_vector_df, dca_proj, flipped))

}


get_classification_accuracy <- function(sample_scores, positive_val){

    if(sum(colnames(sample_scores) %in% c("sample_id", "pathway_score")) != 2){
        stop("sample_scores need column names pathway_score and sample_id")
    }

    pred_dca = ROCR::prediction(sample_scores$pathway_score, sample_scores$Y == positive_val)
    perf_dca_roc = ROCR::performance(pred_dca, "tpr", "fpr")
    auc_dca = ROCR::performance(pred_dca, "auc")
    auc_roc = unlist(auc_dca@y.values)

    perf_dca_pr = ROCR::performance(pred_dca, "prec", "rec")

    x <- perf_dca_pr@x.values[[1]] # Recall values
    y <- perf_dca_pr@y.values[[1]]

    auc_pr = try(MESS::auc(x,y, type = 'spline'), TRUE)
    if(inherits(auc_pr, "try-error")){
        auc_pr = NA
        warning("AUC for PR curve could not be calculated")
    }


    return(list(auc_pr=auc_pr, auc_roc=auc_roc, perf_pr=perf_dca_pr, perf_roc=perf_dca_roc))

}

get_new_samp_score <- function(gene_weights, expression_matr, gene_ids, run_normalization=T){

    if(sum(colnames(expression_matr) %in% c("sample_id")) != 1){
        stop("Need column name sample_id")
    }

    if(sum(colnames(gene_weights) %in% c("gene_weight", "gene_ids")) != 2){
        stop("Need column names gene_weight and gene_ids")
    }
    has_Y = "Y" %in% colnames(expression_matr)

    total_gene_ids = setdiff(colnames(expression_matr), c("Y", "sample_id"))
    gene_ids = intersect(gene_ids, gene_weights$gene_ids)

    if(length(gene_ids) != gene_weights$gene_ids){
        warning("Genes missing in gene_weights or expression_matr")
    }
    # normalize
    # normalize
    dca_matr = expression_matr[,total_gene_ids]
    dca_data = as.matrix(t(scale(log10(t(dca_matr+1)))))
    dca_data[is.nan(dca_data)] <- 0
    row.names(dca_data) = expression_matr$sample_id

    # must preserve ordering of proj_vector
    dca_data = dca_data[,gene_ids]

    # format the projection vector
    proj_vector = gene_weights$gene_weight
    names(proj_vector) = gene_weights$gene_ids
    proj_vector = proj_vector[gene_ids]

    if(sum(names(proj_vector) != colnames(dca_data)) > 0){
        print(names(proj_vector))
        print(colnames(dca_matr))
        stop("Error Matching Gene ids")
    }

    # get projection
    dca_proj = dca_data %*% proj_vector
    dca_proj = data.frame(pathway_score=dca_proj, sample_id=row.names(dca_proj))
    dca_proj = merge(expression_matr, dca_proj, by="sample_id")
    dca_proj = dca_proj[order(dca_proj$pathway_score),]

    if(has_Y){
        dca_proj = dca_proj[,c("sample_id", "Y", "pathway_score")]
    }else{
        dca_proj = dca_proj[,c("sample_id", "pathway_score")]
    }

    return(dca_proj)

}
