\name{get_new_samp_score}
\alias{get_new_samp_score}

\title{
Get a pathway score for an unseen sample
}
\description{
Using the gene weights learned from the reference cohort, we apply the
weightings to new samples to estimate their pathway activity.
}
\usage{
get_new_samp_score(gene_weights, expression_se, gene_ids, run_normalization = TRUE)
}

\arguments{
  \item{gene_weights}{
    This is a data.frame containing gene ids and gene weights, output by
    get_gene_weights. The gene ids must be in the column ids of expression_matr.

}
  \item{expression_se}{
   This is an SummarizedExperiment object of the reference samples. Rows are
   genes and columns are samples. The colData component must contain columns
   \code{Y} and \code{sample_id}. The former indicates whether this is a
   positive or negative sample and the latter is the unique sample id.  
   Within this method, there is a normalization step where each sample is scaled
   across all genes in the SummarizedExperiment assay. For this to be stable and
   consistent, we recommend that the assay contain at least 500 genes that are 
   consistently expressed across all samples in addition to the genes in the 
   pathway of interest.

}
  \item{gene_ids}{
   This is a vector of strings, where each element is a \code{gene_id} in the 
   pathway of interest. The \code{gene_id}s must be present in
   \code{rownames(expression_se)}. 
}

  \item{run_normalization}{
   Boolean value. If TRUE, the data will be log-transformed, centered and scaled.
   This is recommended since this is done to the reference set when learning the
   gene weights.
}
}


\value{
    A data.frame containing the sample id, sample score, and associated Y value
    if it was included in expression_se.
}


\author{
Natalie R. Davidson
}


\examples{
data(tcga_expr_df)

# transform from data.frame to SummarizedExperiment
tcga_se <- SummarizedExperiment(t(tcga_expr_df[ , -(1:4)]),
                                colData=tcga_expr_df[ , 2:4])
colnames(tcga_se) <- tcga_expr_df$tcga_id
colData(tcga_se)$sample_id <- tcga_expr_df$tcga_id

# get the genes of interest, here hypoxia genes
hypoxia_gene_ids <- get_hypoxia_genes()
hypoxia_gene_ids <- intersect(hypoxia_gene_ids, rownames(tcga_se))

# label the samples for classification
colData(tcga_se)$Y <- ifelse(colData(tcga_se)$is_normal, 0, 1)

# now we can get the gene weightings
res <- get_gene_weights(tcga_se, hypoxia_gene_ids, unidirectional=TRUE)
gene_weights <- res[[1]]
sample_scores <- res[[2]]

# get the new data so we can apply our score to it
data(new_samp_df)
new_samp_se <- SummarizedExperiment(t(new_samp_df[ , -(1)]),
                                    colData=new_samp_df[ , 1, drop=FALSE])
colnames(colData(new_samp_se)) <- "sample_id"

new_score_df_calculated <- get_new_samp_score(gene_weights, new_samp_se)

}
