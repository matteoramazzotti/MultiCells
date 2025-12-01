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
- a strategy to query NCBI GEO datasets for specific cells line-based experiments
- a function that retrieves sample names from the selected experiments to allow a proper evaluation of cell identity
- a function to download and assemble RNA-seq data exactly corresponding to the selected samples from multiple GEO series (e.g. subsets of intersting samples only)

It is useful to remind that RNA-seq data form NCBI are not user supplied but are **uniformly remapped, quantified and annotated against a common referece genome**, currently the Human.GRCh38.p13.

After such preparation a data matrix is returned to analyzed with

**groupCapture.R**

an interactive tool to capture groups of samples identified as points in a 2d scatterplot produced by dimensionality reduction techniques such as e.g. PCA, tSNE, UMAP, NMDS , PCoA 

This interactive approach is needed since the groups to be formed and compared are not predetermined (in fact all samples are supposed to be identical), rather they arise upon observation. Basically, the user draws polygons around interesing groups and the script returns data frames containing ready to use data matrices for testing. So easy.

Being RNA-seq count data, the obtained data can be analyzed with popular packages such as DESeq, EdgeR, Limma and so on. **This is up to you**

## Pipeline Example

### 1a - **1.query_performer.pl**

Start by making an index of the available GEOs, we will use a query to fetch all the high throughput sequencing MCF7 related experiments: 

```bash
query='"expression profiling by high throughput sequencing"[Dataset Type] AND "MCF7"[All Fields] AND "rnaseq counts"[Filter] AND "Homo sapiens"[porgn] AND "breast"[All Fields]"'

perl 1.query_performer.pl -q "${query}" -p your_name_for_files_here -O destination_folder
```

The output index will be stored in the *destination_folder/indexes* folder and will have a name like *query_your_name_for_files_here*. Ideally, *destination_folder* should be named after the cell line considered, e.g. "MCF7".

### 1b - **Select Series of Interest**

Select the GSMs to keep (e.g. Control samples) and with them make the *index.tsv.selected* file (a manually curated subset of the previous index that we have downloaded). 
To help yourself in this step, it is handy to do:

```bash
ALIAS="MDAMB231|MDA-MB-231"
grep -i --perl-regexp $ALIAS  MDAMB231/indexes/query_mdamb231.tsv | grep -i --perl-regexp "ctrl|control" > query_mdamb231.subset.tsv
```
> [!NOTE]  
> Edit the ALIAS variable to accomodate keywords compatible with your cell line separated by "|" characters. The `grep` command with the `-i` flag is case insensitive, so don't worry about lower and upper case letters.

Then read the content of this file and select the experiments to keep, storing them into the *index.tsv.selected* file. In the following lines I assume that you store it in the indexes folder.


Make the selected_gse.txt file, a list of unique GSEs from the downladed index:

```bash
awk -F "\t" '{print $2}' destination_folder/indexes/index.tsv.selected | sort | uniq > destination_folder/selected_gse.txt
```

Make the selected_gsm.txt file, a list of unique GSMs from the downloaded index:

```bash
awk -F "\t" '{print $3}' destination_folder/indexes/index.tsv.selected | sort | uniq > destination_folder/selected_gsm.txt
```

Make the gsm_to_gse_selected.tsv file, a list of unique GSMs related to their GSE from the downloaded index:

```bash
awk -F"\t" '{print $3"\t"$2}' destination_folder/indexes/index.tsv.selected > destination_folder/gsm_to_gse_selected.tsv
```

### 2 - **2.get_matrixes.pl**

Download the count matrixes:

```bash
perl 2.get_matrixes.pl -l destination_folder/selected_gsm.txt -L destination_folder/selected_gse.txt -O destination_folder 

```

### 3 - **3.cross_annotation.pl**

Convert entrez gene IDs to gene symbols:

```bash
perl 3.cross_annotation.pl destination_folder/raw_counts/matrix.tsv > destination_folder/matrix_symbol.tsv
```

### 4 - **4.pca_cluster_selector.R**

In R environment, copy istructions from 4.pca_cluster_selector.R. Be sure to set all the paths to the desired folder (i.e. the *destination_folder* from earlier sections).
You will be prompted to select groups in the pca plot. To increase the number of groups to select, edit the *g* argument in the *groupCapture* function.
You will obtain timestamp marked files: e.g. selected_data_Nov_26_2025_09_11_10.tsv and selected_metadata_Nov_26_2025_09_11_10.tsv files.


### 5 - **5.de.R**


Perform differential expression analysis using timestamp marked files obtain previously:

```bash
Rscript 5.de.R -c MCF7/selected_data_Nov_26_2025_09_11_10.tsv -m MCF7/selected_metadata_Nov_26_2025_09_11_10.tsv -M MCF7
```

You will obtain more timestamp marked files: e.g. de_data_Nov_26_2025_09_11_10.tsv and de_data_subset_Nov_26_2025_09_11_10.tsv files.


### 6 - **6.volcano_plot.R**

Obtain volcano plot based on total de data:

```bash
Rscript 6.volcano_plot.R -i MCF7 -R MCF7/de_data_Nov_26_2025_09_11_10.tsv -T 1.5 -t 0.05
```

`-T 1.5` stands for a log2FC threshold of 1.5 and `-t 0.05` stands for a padj threshold of 0.05. Plots are stored in the plots folder.


### 7 - **7.functional_analysis.R**

Finally, we perform an enrichment analysis using clusterProfiler:

```bash
Rscript clusterProfiler.R --genes MCF7/de_data_subset_Nov_26_2025_09_11_10.tsv --gmt gmt/c5.go.bp.v2025.1.Hs.symbols.gmt --out-prefix c5.go.bp --out-dir MCF7
```

Be sure to use the `de_data_subset_` file with the `--genes` argument. The `--out-prefix` argument is used to customize the name of the output file, we are using one based on the gmt employed for this analysis, so the final output name will be `cp_results_c5.go.bp.tsv`. This will also return a barplot with the top 20 pathways enriched.
