\name{get_hypoxia_genes}
\alias{get_hypoxia_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Ensembl ids of hypoxia related genes.
}
\description{
Returns a vector of Ensembl ids of hypoxia related genes.
}
\usage{
get_hypoxia_genes()
}

\value{
Vector of ensembl ids.
}

\author{
Natalie R. Davidson
}


\examples{


# read in the reference expression data for hypoxia score generation
data(tcga_expr_df)

# transform from data.frame to SummarizedExperiment
tcga_se <- SummarizedExperiment(t(tcga_expr_df[ , -(1:4)]),
                                colData=tcga_expr_df[ , 2:4])
colnames(tcga_se) <- tcga_expr_df$tcga_id
colData(tcga_se)$sample_id <- tcga_expr_df$tcga_id

# let's get the expression of hypoxia associated genes
hypoxia_gene_ids <- get_hypoxia_genes()
hypoxia_gene_ids <- intersect(hypoxia_gene_ids, rownames(tcga_se))
hypoxia_se <- tcga_se[hypoxia_gene_ids,]


}
