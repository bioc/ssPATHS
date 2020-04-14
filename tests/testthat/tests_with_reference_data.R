test_that("Gene weight on reference dataset works", {

    data(tcga_expr_df)

    # transform from data.frame to SummarizedExperiment
    tcga_se <- SummarizedExperiment(t(tcga_expr_df[ , -(1:4)]),
                                    colData=tcga_expr_df[ , 2:4])
    colnames(tcga_se) <- tcga_expr_df$tcga_id
    colData(tcga_se)$sample_id <- tcga_expr_df$tcga_id


    hypoxia_gene_ids <- get_hypoxia_genes()
    hypoxia_gene_ids <- intersect(hypoxia_gene_ids, rownames(tcga_se))
    hypoxia_se <- tcga_se[hypoxia_gene_ids,]

    colData(hypoxia_se)$Y <- ifelse(colData(hypoxia_se)$is_normal, 0, 1)

    # now we can get the gene weightings
    res <- get_gene_weights(hypoxia_se)
    gene_weights_test <- res[[1]]
    sample_scores <- res[[2]]

    data("gene_weights_reference")
    expect_equal(unlist(gene_weights_test[,1]), unlist(gene_weights_reference[,1]))
    expect_equal(gene_weights_test[,2], as.character(gene_weights_reference[,2]))
    expect_equal(row.names(gene_weights_test), row.names(gene_weights_reference))
})

test_that("test classification on reference dataset works", {

    data(tcga_expr_df)

    # transform from data.frame to SummarizedExperiment
    tcga_se <- SummarizedExperiment(t(tcga_expr_df[ , -(1:4)]),
                                    colData=tcga_expr_df[ , 2:4])
    colnames(tcga_se) <- tcga_expr_df$tcga_id
    colData(tcga_se)$sample_id <- tcga_expr_df$tcga_id


    hypoxia_gene_ids <- get_hypoxia_genes()
    hypoxia_gene_ids <- intersect(hypoxia_gene_ids, rownames(tcga_se))
    hypoxia_se <- tcga_se[hypoxia_gene_ids,]

    colData(hypoxia_se)$Y <- ifelse(colData(hypoxia_se)$is_normal, 0, 1)

    # now we can get the gene weightings
    res <- get_gene_weights(hypoxia_se)
    sample_scores <- res[[2]]

    training_res <- get_classification_accuracy(sample_scores, positive_val=1)
    print(training_res[[2]]-0.91112)

    expect_equal(round(training_res[[2]], 5), 0.91112)
})


test_that("test classification on new dataset works", {

    data(tcga_expr_df)

    # transform from data.frame to SummarizedExperiment
    tcga_se <- SummarizedExperiment(t(tcga_expr_df[ , -(1:4)]),
                                    colData=tcga_expr_df[ , 2:4])
    colnames(tcga_se) <- tcga_expr_df$tcga_id
    colData(tcga_se)$sample_id <- tcga_expr_df$tcga_id

    hypoxia_gene_ids <- get_hypoxia_genes()
    hypoxia_gene_ids <- intersect(hypoxia_gene_ids, rownames(tcga_se))
    hypoxia_se <- tcga_se[hypoxia_gene_ids,]

    colData(hypoxia_se)$Y <- ifelse(colData(hypoxia_se)$is_normal, 0, 1)

    # now we can get the gene weightings
    res <- get_gene_weights(hypoxia_se)
    gene_weights <- res[[1]]
    sample_scores <- res[[2]]

    data(new_samp_df)
    new_samp_se <- SummarizedExperiment(t(new_samp_df[ , -(1)]),
                                        colData=new_samp_df[, 1, drop=FALSE])
    colnames(colData(new_samp_se)) <- "sample_id"

    new_score_df_calculated <- get_new_samp_score(gene_weights, new_samp_se)

    data(expected_score_output)
    row.names(expected_score_output) <- NULL
    new_score_df_calculated <- as.data.frame(new_score_df_calculated)
    row.names(new_score_df_calculated) <- NULL
    expect_equal(new_score_df_calculated, expected_score_output)

})

