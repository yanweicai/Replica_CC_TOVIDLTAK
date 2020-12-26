This is a code repository for the manuscript 'Replicability and reproducibility of genetic analysis between different mice liver studies using the Collaborative Cross'

# Study Introduction
Replicability of experiments is an essential part of the scientific method. By combining research records from several large scale genetics experiments in panels of inbred mice strains (the Collabrative Cross, CC) with mostly identical designs, we have a dataset that to a large extent mimics the process of strict replication.
From this dataset, we applied replicability comparison of both phenotypic and expression data.

![Screenshot](introfigure.png)
Three studies were originally designed to invest drug induced liver injury of three different drugs (Coded as: TOV, TAK and IDL study). 
The control CC mice from the three studies underwent very similar experiments of vehicle-only dosing. The uniformity of experiments, together with genetically identical inbred subjects, makes these studies been considered independent replicating experiments on identical subjects.
Body weight and liver related phenotype are measured. Microarray expression is measured in the TOV and TAK data.

# Files
* phenotype.R
 * Analysis of variance components of phenotypes (Table 1)
 * QTL mapping analysis of phenotypes (Including Figure 2)
 * Analysis of within-strain differences (Figure 6)
* ExpAnalysis.R
 * Variance components analysis in expression (Figure 3)
 * Analysis and plotting for Figure 4,5
* Cluster_ExpMapping_cluster.r
 * eQTL mapping of expression data (Figure 5)
* Cluster_ExpggLASSO.r
 * Lasso, gLasso, ridge, elastic-net analysis of expression data (Figure 4D)
