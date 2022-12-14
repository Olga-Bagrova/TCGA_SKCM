library(tidyverse)
library(gtsummary)

ck.cc <- clusterProfiler::compareCluster(geneCluster = list("1500" = colnames(mat)), 
                                         fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "CC",
                                         pvalueCutoff = 0.01, qvalueCutoff =  0.05)

ck.mf <- clusterProfiler::compareCluster(geneCluster = list("1500" = colnames(mat)), 
                                         fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "MF",
                                         pvalueCutoff = 0.01, qvalueCutoff =  0.05)

ck.bp <- clusterProfiler::compareCluster(geneCluster = list("1500 genes" = colnames(mat)),
                                         fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP",
                                         pvalueCutoff = 0.01, qvalueCutoff =  0.05)


clusterProfiler::dotplot(ck.bp, by = "count", showCategory = 20, 
                         title = "Gene Ontology enrichment") + 
  scale_colour_gradient(limits=c(0, 0.001), low="red", high="blue") +
  theme(plot.title = element_text(hjust = 0.5))

