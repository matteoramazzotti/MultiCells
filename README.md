# MultiCells - "Single Bulk" Cell Analysis in Linux and R

MultiCells is a project built upon the necessity of **evaluating the real differences in cell lines supposed to be almost identical**.

We devise this **"single bulk" cell analysis** as a way to compare the RNA-seq based expression profiles among cell lines from different experiments.  

**Possible uses**
- Check consistency across untrated cell lines from different experiments 
- Establish a robust baseline experiment
- Assess the effect of a drug or a treatment in multiple cells from different experiments
- ...

The **MultiCells** approach exploits the NCBI GEO database to download bulk RNA-seq data from cultured cell lines.
It contains:
- a strategy to query NCI GEO datasets for specific cells line-based experiments
- a function that retrieves sample names from the selected experients to allow a proper evaluation of cell identity
- a function to download and assemble RNA-seq data exactly corresponding to the selected samples from multiple GEO series ((e.g. subsets of intersting samples only)

It is useful to remind that RNA-seq data form NCBI are not user supplied but are **uniformly remapped, quantified and annotated against a common referece genome**, currently the Human.GRCh38.p13.

After such preparation a data matrix is returned to analyzed with

**groupCapture.R**

an interactive tool to capture groups of samples identified as points in a 2d scatterplot produced by dimensionality reduction techniques such as e.g. PCA, tSNE or nmds 

This interactive approach is needed since the groups to be formed and compared are not predetermined (e.g. all samples are supposed to be identical), 
rather they arise upon observation. Basically, the user draws polygons around interesing groups and the script returns data frames containing 
ready to used data matrices for testing. So easy.

Being RNA-seq count data, the obtained data can be analyzed with popular packages such as DESeq, EdgeR, Limma and so on. **This is up to you**
