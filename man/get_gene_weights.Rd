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
get_gene_weights(expression_se, gene_ids, unidirectional)
}

\arguments{
  \item{expression_se}{
   This is an SummarizedExperiment object of the reference samples. Rows are
   genes and columns are samples. The colData component must contain a 
   \code{sample_id} column. Within this method, there is a normalization step 
   where each sample is scaled across all genes in the SummarizedExperiment 
   assay. For this to be stable and consistent, we recommend that the assay 
   contain at least 500 genes that are consistently expressed across all samples
   in addition to the genes in the pathway of interest.
}
  \item{gene_ids}{
   This is a vector of strings, where each element is a \code{gene_id} in the 
   pathway of interest. The \code{gene_id}s must be present in
   \code{rownames(expression_se)}. 
}
  \item{unidirectional}{
   This is a boolean, \code{default=TRUE}. Most genesets are unidirectional, 
   meaning that most genes are either increasing or decreasing together. If this
   is set to  \code{TRUE}, then the learned weights will be clipped such that 
   the dominant directionality is kept, and the other gene weights are set to 
   zero.
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

# setup labels for classification
colData(tcga_se)$Y <- ifelse(colData(tcga_se)$is_normal, 0, 1)

# now we can get the gene weightings
res <- get_gene_weights(tcga_se, hypoxia_gene_ids)
gene_weights_test <- res[[1]]
sample_scores <- res[[2]]


}
