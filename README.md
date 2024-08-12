# RiskGenesForCKD
Scripts for single-cell analysis part for manuscript "Multi-Scalar Data Integration Decoding Risk Genes for Chronic Kidney Disease"

1. preprocess.R
   The preprocess of data downloaded from open source provided by authors.
   Raw counts and original sample information were used. 
   Normalisation processes by standard Seurat pipeline.
   Orginal annotations were retained and adopted in the Seruat object. 

   Data set 1. publication  https://zenodo.org/record/4059315
   Data set 2. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211785

2. integrate. R
   Integration of two kidney data sets that were preprocessed by 1.preprocess.R
   Raw plots for the manuscript.  
