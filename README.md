# MultiCells - "Single Bulk" Cell Analysis in Linux and R

MultiCells is a project built upon the necessity of **evaluating the real differences in cell lines supposed to be almost identical**.

We devise this **"single bulk" cell analysis** as a way to compare the RNA-seq based expression profiles among cell lines from different experiments.  

## **Possible uses**
- Check consistency across untrated cell lines from different experiments 
- Establish a robust baseline experiment
- Assess the effect of a drug or a treatment in multiple cells from different experiments
- ...

The **MultiCells** approach exploits the NCBI GEO database to download bulk RNA-seq data from cultured cell lines.
It contains:
- a strategy to query NCI GEO datasets for specific cells line-based experiments
- a function that retrieves sample names from the selected experients to allow a proper evaluation of cell identity
- a function to download and assemble RNA-seq data exactly corresponding to the selected samples from multiple GEO series (e.g. subsets of intersting samples only)

It is useful to remind that RNA-seq data form NCBI are not user supplied but are **uniformly remapped, quantified and annotated against a common referece genome**, currently the Human.GRCh38.p13.

After such preparation a data matrix is returned to analyzed with

**groupCapture.R**

an interactive tool to capture groups of samples identified as points in a 2d scatterplot produced by dimensionality reduction techniques such as e.g. PCA, tSNE or nmds 

This interactive approach is needed since the groups to be formed and compared are not predetermined (in fact all samples are supposed to be identical), 
rather they arise upon observation. Basically, the user draws polygons around interesing groups and the script returns data frames containing 
ready to use data matrices for testing. So easy.

Being RNA-seq count data, the obtained data can be analyzed with popular packages such as DESeq, EdgeR, Limma and so on. **This is up to you**

## Pipeline Example

### 1a - **1.query_performer.pl**

Start by making an index of the available GEOs, we will use a query to fetch all the high throughput sequencing MCF7 related experiments: 

```
	query='"expression profiling by high throughput sequencing"[Dataset Type] AND "MCF7"[All Fields] AND "rnaseq counts"[Filter] AND "Homo sapiens"[porgn] AND "breast"[All Fields]"'

	perl 1.query_performer.pl -q "${query}" -p your_name_for_files_here
```

The output index will be stored in the indexes folder and will have a name like query_your_name_for_files_here

### 1b - **Select Series of Interest**

Select the GSMs to keep (e.g. Control samples) and with them make the index.tsv.selected file (a manually curated subset of the previous index that we have downloaded). In the following lines I assume that you store it in the indexes folder.

Make the selected_gse.txt file, a list of unique GSEs from the downladed index:

```
	awk -F "\t" '{print $2}' indexes/index.tsv.selected | sort | uniq > selected_gse.txt
```

Make the selected_gsm.txt file, a list of unique GSMs from the downloaded index:

```
	awk -F "\t" '{print $3}' indexes/index.tsv.selected | sort | uniq > selected_gsm.txt
```

Make the gsm_to_gse_selected.tsv file, a list of unique GSMs related to their GSE from the downloaded index:

```
awk -F"\t" '{print $3"\t"$2}' indexes/index.tsv.selected > gsm_to_gse_selected.tsv
```

### 2 - **2.get_matrixes.pl**

Download the count matrixes:

```
	perl 2.get_matrixes.pl

```

### 3 - **3.cross_annotation.pl**

Convert entrez gene IDs to gene symbols:

```
	perl 3.cross_annotation.pl raw_counts/matrix.tsv > matrix_symbol.tsv
```

### 4 - **4.pca_cluster_selector.R**

In R environment, copy istructions from 4.pca_cluster_selector.R.
You will be prompted to select groups in the pca plot. To increase the number of groups to select, edit the *g* argument in the *groupCapture* function.
You will obtain selected_data.tsv and selected_data.tsv files.


### 5 - **5.de.R**

In R environment, copy istructions from 5.de.R and perform the differential expression analyses, remember to edit the instructions to accomodate the groups you selected.
