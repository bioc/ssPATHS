\name{get_classification_accuracy}
\alias{get_classification_accuracy}

\title{
Get Classification Accuracy
}
\description{
Get the AUC-ROC, AUC-PR, and ROC/PR curves for plotting.
}
\usage{
get_classification_accuracy(sample_scores, positive_val)
}

\arguments{
  \item{sample_scores}{
   This is a data.frame containing the sample id, score, and true label Y. This
   object is returned by the method get_gene_weights.
}
  \item{positive_val}{
   This is the value that will denote a true positive. It must be one of the two
   values in the Y column in sample_scores.
}
}

\value{
This returns a list of performance metrics
  \item{auc_pr }{Area under the PR-curve}
  \item{auc_roc }{Area under the ROC-curve}
  \item{perf_pr }{ROCR object for plotting the PR-curve}
  \item{perf_roc }{ROCR object for plotting the ROC-curve}

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

hypoxia_gene_ids <- get_hypoxia_genes()
hypoxia_gene_ids <- intersect(hypoxia_gene_ids, rownames(tcga_se))

colData(tcga_se)$Y <- ifelse(colData(tcga_se)$is_normal, 0, 1)

# now we can get the gene weightings
res <- get_gene_weights(tcga_se, hypoxia_gene_ids, unidirectional=TRUE)
sample_scores <- res[[2]]

# check how well we did
training_res <- get_classification_accuracy(sample_scores, positive_val=1)
print(training_res[[2]])

plot(training_res[[3]], col="orange", ylim=c(0, 1))
legend(0.1,0.8,c(training_res$auc_pr,"\n"), border="white", cex=1.7,
        box.col = "white")

plot(training_res[[4]], col="blue", ylim=c(0, 1))
legend(0.1,0.8,c(training_res$auc_roc,"\n"),border="white",cex=1.7,
        box.col = "white")

}
