\name{get_new_samp_score}
\alias{get_new_samp_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get a pathway score for an unseen sample
}
\description{
Using the gene weights learned from the reference cohort, we apply the
weightings to new samples to estimate their pathway activity.
}
\usage{
get_new_samp_score(gene_weights, expression_matr, run_normalization = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_weights}{
    This is a data.frame containing gene ids and gene weights, output by
    get_gene_weights. The gene ids must be in the column ids of expression_matr.

}
  \item{expression_matr}{
   This is an expression matrix of 1 or more new samples. Rows are
   samples and mandatory columns are sample_id (unique sample id) and gene
   expression. If a "Y" column is included, it is ignored in the calculation but
   included in the output. All columns not titled "Y" or "sample_id" are assumed
   to be gene names in the pathway of interest.

}
  \item{run_normalization}{
   Boolean value. If True, the data will be log-transformed, centered and scaled.
   This is recommended since this is done to the reference set when learning the
   gene weights.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    A data.frame containing the sample id, sample score, and associated Y value
    if it was included in expression_matr.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Natalie R. Davidson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library("dml")
library("data.table")
library("ROCR")
library("MESS")

# read in the reference expression data for hypoxia score generation
data(tcga_expr_df)

# take a look at the format, we need a Y column and a sample_id column.
tcga_expr_df[1:10,1:10]

# let's get the info and columns we need to learn a hypoxia score
hypoxia_gene_ids = get_hypoxia_genes()
hypoxia_gene_ids = intersect(hypoxia_gene_ids, colnames(tcga_expr_df))
hypoxia_df = tcga_expr_df[,c("tcga_id", "is_normal", hypoxia_gene_ids)]
colnames(hypoxia_df)[1:2] = c("sample_id", "Y")

# format the Y column such that normal samples are 0 and tumor samples are 1
hypoxia_df$Y = 0
hypoxia_df$Y[tcga_expr_df$is_normal==TRUE] = 0
hypoxia_df$Y[tcga_expr_df$is_normal==FALSE] = 1

# now we can get the gene weightings
res = get_gene_weights(hypoxia_df)
gene_weights = res[[1]]
sample_scores = res[[2]]

# Now, using the gene weightings learned from the reference set, we can apply
# it to a new sample.
data("new_samp_df")
new_score_df = get_new_samp_score(gene_weights, new_samp_df)
new_score_df

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
