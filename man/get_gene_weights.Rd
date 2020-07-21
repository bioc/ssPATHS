\name{get_gene_weights}
\alias{get_gene_weights}

\title{
  Get Gene Weights from Reference Data
}
\description{
  This method performs linear discriminant analysis on a reference dataset
  using a pre-defined set of genes related to a pathway of interest.
}
\usage{
get_gene_weights(expression_se)
}

\arguments{
  \item{expression_se}{
   This is an SummarizedExperiment object of the reference samples. Rows are
   genes and columns are samples. The colData component must contain columns
   \code{Y} and \code{sample_id}. The former indicates whether this is a
   positive or negative sample and the latter is the unique sample id.
   All genes in the SummarizedExperiment object are assumed to be genes
   in the pathway of interest.
}
}

\value{
  A list containing the gene weights and estimated scores of the reference
  samples.
  \item{proj_vector_df }{A dataframe containing the gene weights and gene ids}
  \item{dca_proj }{A dataframe containing the sample scores and sample ids.}
}
\references{
 Steven C.H. Hoi, W. Liu, M.R. Lyu and W.Y. Ma (2006). Learning Distance
 Metrics with Contextual Constraints for Image Retrieval. Proceedings IEEE
 Conference on Computer Vision and Pattern Recognition (CVPR2006).
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

# get related genes, for us hypoxia
hypoxia_gene_ids <- get_hypoxia_genes()
hypoxia_gene_ids <- intersect(hypoxia_gene_ids, rownames(tcga_se))
hypoxia_se <- tcga_se[hypoxia_gene_ids,]

# setup labels for classification
colData(hypoxia_se)$Y <- ifelse(colData(hypoxia_se)$is_normal, 0, 1)

# now we can get the gene weightings
res <- get_gene_weights(hypoxia_se)
gene_weights_test <- res[[1]]
sample_scores <- res[[2]]


}
