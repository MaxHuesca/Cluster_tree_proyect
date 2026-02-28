# Creating dendrogram trees using blast alignments by diferent hierarchical cluster methods  
>Clustering module Bioinforamtics II
>
>Ismael Maximiliano De Los Santos Huesca

## Summary 

Hierarchical clustering methods offer a framework for exploring relationships within biological datasets by generating dendrograms that represent distance-based groupings. When applied to protein sequence data, these approaches can reveal structural and functional similarities without requiring prior phylogenetic information.

We performed all‑vs‑all BLASTP alignments on three protein datasets: 100 ABC family proteins (expected 4 clusters), 100 H1 histones from four taxonomic groups (*Mammalia*, *Protostomia*, *Fungi*, *Echinodermata*), and 100 cytochrome C sequences with no prior grouping information. Alignment scores were normalized and converted to dissimilarity matrices (1 − normalized score). Cluster numbers were inferred using Elbow (WSS), Silhouette, and Gap statistic methods across four merging criteria (ward.D2, single, complete, average). Dendrograms were constructed using agglomerative (hclust, agnes) and divisive (diana) algorithms and compared for topological informativeness.

The `Silhouette` method provided the most consistent cluster number inferences across all datasets, while Gap statistic frequently yielded inconclusive results. Ward merging consistently produced the most interpretable dendrograms with clearly separated branches, whereas single linkage generated poorly resolved trees with the lowest agglomerative coefficients. For ABC proteins, ward‑based dendrograms recovered the expected four clusters; for H1 histones, they suggested 2‑3 groups, potentially reflecting closer affinity between Protostomia and Echinodermata sequences. Cytochrome C analysis indicated two major clusters, consistent with Silhouette inference. Notably, dendrograms with the highest agglomerative coefficients (e.g., single linkage) were often the least informative biologically.

## Contents 
This proyec focus in the use of **hierarchical cluster methods** in order to create a dendrogram, this is relevany beacause we are using as input data the score results froma a all vs all 
blastp alignments, this analysis is structure in the following way: 
* [`Data`](https://github.com/MaxHuesca/Cluster_tree_proyect/tree/main/data): Directory with all the raw and clean data, it has the fasta sequences and the alignments results
* [`results`](https://github.com/MaxHuesca/Cluster_tree_proyect/tree/main/results): Here are all the final results for the scripts that are in source, at the firts level the .tsv and matrix results are found
    * [`plots`](https://github.com/MaxHuesca/Cluster_tree_proyect/tree/main/results/plots): Here are stored all the images that the analysis generate, at the firts level you can find the inference methods plots 
        * [`dendrograms`](https://github.com/MaxHuesca/Cluster_tree_proyect/tree/main/results/plots/dendrograms): The dendrograms are the final result of the analyisis an here they are stored
    * [`trees`](https://github.com/MaxHuesca/Cluster_tree_proyect/tree/main/results/trees): Here are stored all the newick trees 
* [`docs`](https://github.com/MaxHuesca/Cluster_tree_proyect/tree/main/docs): This sirectory are for format files and also binaries ones 
* [`src`](https://github.com/MaxHuesca/Cluster_tree_proyect/tree/main/src): Here are all the steps of the analysis in diferent scripts that converge in the [`proyect_report.qmd`](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/src/proyect_report.qmd) in the following order of execution:
    1. [`create_tabular.sh`](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/src/create_tabular.sh): This script has as input the fasta file with all the sequences and the number of sequences the user wants to retrive and thus analyze, it generates a .tsv file with the results.
    2. [`parsed_tabular.py`](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/src/parsed_parseTabular.py): This script make the disimilarity matrix and in the process generates the score matrix and the normalized matrix, the input is the previous generate .tsv file an the following options to make the matrix:
    ```bash 
    options:
  -h, --help       show this help message and exit
  -t, --tsv TSV    Path of the tsv file with the results of the alignment in the format query /t subject /t score
  -w, --way WAY    Specify the way to build the matrix 'upp' for consider the query result over the subject (upper
                   triangular), 'max' for consider the max result
  -o, --outB OUTB  baseline for the out format
    ```
    3. [`proyect_report.qmd`](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/src/proyect_report.qmd): This qmd has all the summary of the process and make the formal clustering analysis, the inputs are the disimilarity matrixes.
    4. [`cluster_tree.qmd`](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/src/cluster_tree.qmd):This qmd has onky the tree proccesing usign no ggplot objects to make the trees. 

```bash
├── data/               # Raw and cleaned sequences, BLAST results
├── results/            
│   ├── plots/          # Inference plots and dendrograms
│   └── trees/          # Newick format trees
├── src/                
│   ├── create_tabular.sh           # BLAST execution
│   ├── parsed_parseTabular.py       # Matrix construction
│   ├── proyect_report.qmd           # Main analysis
│   └── cluster_tree.qmd             # Tree-only processing
└── docs/               # Documentation files
```

## Introduction

With the creation of human knowledge, it came out a necessity for clustering things and objects in order to understand in a better way the world and things that human beings see. That process involves selecting attributes from the objects that can discriminate them from others, so for that, it is very important to select those characteristics that can differentiate the parts of a set. Thus, that is one of the major challenges when we talk about the task of clustering things.

Biology does not go far dealing with that problem; branches of this discipline face this problem every day. Things like the taxonomy features need to be precisely assigned in order to get into different sets all the living beings that biology can recover. Actually, something that has an enormous deal with that is all the biology tasks related to novel techniques like NGS applied to molecules like DNA and RNA, and also branches of the same tree related to other types of biomolecules like proteins (Proteomics).

The objective of this project is to make a dendrogram by applying hierarchical clustering methods. This procedure allows us to do so because, different from the non-hierarchical algorithms, it solves all the possible number of cluster solutions for a given data set, same that can be related to its values with distances from each other. So these solutions provided can be represented as dendrograms where the branches are the distances between all the pairs. 

---

## Results  
### Selecting and parsing the data

For this analysis we are going to use 3 different sets of data because of their nature:
* First, we are analyzing the data set provided by ***Dr. Luis Arturo Medrano Soto*** that contains a set of 100 proteins that belong to the **ABC** family from which we can cluster 4 or 3 groups related to the domains embedded in those proteins.
* In the second group, we have a set of 100 **H1** proteins from 4 different life clades (25 for each one):
    * *Mammalia*
    * *Protostomia*
    * *Fungi*
    * *Echinodermata*
* Finally, we have a set of 100 **cytochrome C** proteins for which we do not have prior phylogenetic information.

The two groups that were not provided by the course titular were extracted from [UniProt](https://www.uniprot.org/) and [EBI](https://www.ebi.ac.uk/interpro/) databases respectively. They were parsed for cleaning their headers using `awk` from the command line.
```{bash}
#| eval: false

awk 'BEGIN{FS="|"} /^>/{print ">" $2; next} {print}' ../data/H1_seqs_raw.faa > ../data/H1_seqs.faa

awk '/^>/{sub(/\|.*/,"",$0)}1' ../data/citrocromeC_seq_hmmr.faa > ../data/citrocromeC_seq_hmmr_clean.faa
```

After cleaning the headers, they are ready for the first step. 

### Blasting

To access "distances" in our data, we perform alignments in an all-vs-all format. This was done with the suite [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) using output format 7 in order to retain only the query, subject, and score columns. This is accessed through the script `src/create_tabular`; for more information, check this script. With this data, we can parse it to generate a matrix. The table can be seen as follows:

| query   | subject | score |
|---------|---------|-------|
| Q27749  | P02290  | 258.0 |
| Q27749  | P06145  | 195.0 |
| Q27749  | P06146  | 188.0 |
| Q27749  | P16887  | 187.0 |
| Q27749  | P07794  | 181.0 | 

Once we have this TSV file, we can start with the matrix construction.

---

###  Matrix Construction

The task of building the matrix is handled by the script `src/parsed_parseTabular.py`, which ensures that the matrix is constructed in a symmetric way and in two possible ways. First, considering the alignment where protein A is the query over the alignment where the same protein is the subject (option `upp` in the flag `-w`), and second, considering the alignment with the maximum score (option `max` in the flag `-w`).

This script also normalizes the matrix using as the normalization factor the maximum score of the matrix without considering the diagonal values (because they correspond to self-self alignments, where the score is maximum and thus the distance is 0). After normalization, it generates a dissimilarity matrix by subtracting the normalized score from 1.
$$
1 - B_{i,j}
$$ 

| query   | Q27749   | P02287   | P02290   | P04735   | P06145   | P06146   | P07794   | P13281   | P16887   | ... | Q8CJI4 | Q8IZA3 | Q8N1T3 | Q8VIK3   | Q91YA2 | Q92522 | Q93079   | Q9ULR3   | P06348 | P40286 |
|---------|----------|----------|----------|----------|----------|----------|----------|----------|----------|-----|--------|--------|--------|----------|--------|--------|----------|----------|--------|--------|
| Q27749  | 1.000000 | 0.788484 | 0.696827 | 0.000000 | 0.770858 | 0.779083 | 0.787309 | 0.938425 | 0.780259 | ... | 0.0    | 0.0    | 0.0    | 0.000000 | 0.0    | 0.0    | 0.793184 | 0.977791 | 0.0    | 0.0    |
| P02287  | 0.788484 | 1.000000 | 0.797885 | 0.969683 | 0.796710 | 0.801410 | 0.779083 | 0.000000 | 0.799060 | ... | 0.0    | 0.0    | 0.0    | 0.000000 | 0.0    | 0.0    | 0.774383 | 0.000000 | 0.0    | 0.0    |
| P02290  | 0.696827 | 0.797885 | 1.000000 | 0.973325 | 0.767333 | 0.775558 | 0.786134 | 0.950176 | 0.777908 | ... | 0.0    | 0.0    | 0.0    | 0.000000 | 0.0    | 0.0    | 0.790834 | 0.978731 | 0.0    | 0.0    |
| P04735  | 0.000000 | 0.969683 | 0.973325 | 1.000000 | 0.975088 | 0.000000 | 0.970623 | 0.000000 | 0.000000 | ... | 0.0    | 0.0    | 0.0    | 0.979671 | 0.0    | 0.0    | 0.970153 | 0.000000 | 0.0    | 0.0    |
| P06145  | 0.770858 | 0.796710 | 0.767333 | 0.975088 | 1.000000 | 0.777908 | 0.782609 | 0.947944 | 0.777908 | ... | 0.0    | 0.0    | 0.0    | 0.000000 | 0.0    | 0.0    | 0.787309 | 0.000000 | 0.0    | 0.0    |

Having the dissimilarity matrix is crucial because it is the milestone for reading these objects with the `as.dist()` function in R, the language in which the following steps are performed.

---

### Clustering

When we talk about clustering, we can mention two main methods: hierarchical clustering and non-hierarchical clustering. In this analysis, we are using the first one to obtain all the clusters in which our data can be classified as a function of our selected variable (the distance in the form of a normalized score from an alignment). However, before performing clustering, the question of how many clusters our dataset can contain always arises. To answer this question, we can try to infer the number of clusters that our data heterogeneity could explain using different methods:

#### Elbow (WSS)

The Elbow method focuses on considering the sum of squared errors (SSE). By plotting these errors present in our data, we can try to identify abrupt changes in the tendency while assuming a number of clusters. It can also be modeled using a hierarchical clustering method, so the abrupt change can be viewed as an elbow in a decreasing plot; that is why the method is named this way. This elbow represents the inferred number of clusters our data contains.

#### Silhouette

This heuristic method consists of evaluating how similar an observation is to its own cluster. This can be accessed using distances and is compared with the distance to other clusters. From this, it obtains a silhouette coefficient. The plot explains the trend, and the peaks in the plot represent the inferred number of clusters. As with the Elbow method, it can use a hierarchical clustering method to model cluster formation.

#### Gap (gap_stat)

This heuristic method, unlike the two previous methods, quantifies the distribution of the variation across clusters by comparing it with the theoretical values of a uniform distribution. This can be done using a bootstrapping method.  

#### Inference in H1
![Inference_H1data](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/results/plots/Inf_multiplot_H1.png)
#### H1

In this set of proteins, we also know that the number of clusters is 4. We can note that the gap statistic method also fails to determine a conclusive number of inferred clusters, because we do not observe a dominant peak across those plots. This differs from the plots relative to the silhouette method, in which we can see more consistent peaks across the values 2 and 3, with the number of 2 clusters showing more evidence in the hierarchical clustering using the single merge method. Finally, we can identify a putative value between 2 and 3 clusters in the Elbow method.

#### Inference in ABC 
![Inference_ABCdata](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/results/plots/Inf_multiplot_ABC.png) 

#### ABC

In this case, we previously knew the number of clusters that this dataset is supposed to have, and this is what we observed in the plots. Note that the first row, corresponding to the **gap** method, does not show conclusive results, because in none of the plots we can see the peak we are searching for. This changes in the **silhouette** plots, where we can see that the number of clusters around 4 is more homogeneous (merge methods *ward.D2*, *complete*). This is also reflected in the elbows observed for the **wss** method (*ward.D2*, *complete*). This gives us a more confident assumption that the ABC dataset has 4 clusters representing those protein sequences.

#### Inference in CitochromeC
![Inference_CitochormeC](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/results/plots/Inf_multiplot_CitochromeC.png)  
#### CytochromeC

This could be the most interesting case because we do not have an a priori number of clusters. As in the two previous groups, we can see that the gap statistic method fails to provide conclusive numbers of clusters, as well as the wss method, where we cannot identify a remarkable elbow across the four plots. However, interestingly, in the silhouette method we can observe a consensus peak at 2 clusters, giving us insights into what we could see in the following dendrograms.


With these plots, we can infer a putative number of clusters in our data. The interpretation depends on the method plotted, and sometimes the method used to merge the clusters (criteria) also plays a role in the shape of the plot.

--- 

### Dendrograms

To construct the dendrograms, we are using a hierarchical clustering method with the objective of representing the distances stored in our distance matrix as a solution of all the possible clusters that these distances represent, as a function of an agglomerative or divisive method and its corresponding merge method.

For this task, we are solving all the possible clusters for the three datasets using two agglomerative methods (`hclust`, `agnes`) and one divisive method (`diana`), together with the previous four merging methods.

#### Dendrograms clusters in H1
![Inference_H1data](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/results/plots/dendrograms/H1dendrograms.png)

#### H1

In the H1 proteins, the most informative dendrogram is also the one that uses the `ward` merge method. In this case, we can discriminate between three major branches, suggesting that there are three groups in our data based on the distances present in our matrix, consistent with the results from the silhouette inference method. These results could be explained by the close phylogenetic distance between the *Protostomia* clade and the *Echinodermata* phylum, but this clustering evidence is not conclusive for such assumptions.

Regarding the least informative dendrograms, there is a clear result: those produced by the `single` merge method. This is also reflected in its non-defined agglomerative coefficient, in contrast with the 0.9592296 coefficient corresponding to the `ward` merge method.

#### Dendrograms clusters in ABC 
![Inference_ABCdata](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/results/plots/dendrograms/ABCdendrograms.png)
#### ABC

We can see that the most informative dendrogram is the one corresponding to the `hclust` method with the `ward` merge strategy. This is because we can better discriminate inter-cluster distances, as the branches are more separated. This is consistent with the number of clusters we were expecting, probably with the closest ones being those corresponding to ABC2 "a" and "b". Therefore, we can interpret this dendrogram as the one that best represents the taxonomy of the proteins.

In contrast, the least informative trees are those corresponding to the `single` merge method, which are also the ones with the lowest agglomerative coefficient (0.4836667). Note that their branches do not diverge enough to identify our initial four groups.

#### Dendrograms clusters in CitochromeC
![Inference_CitochormeC](https://github.com/MaxHuesca/Cluster_tree_proyect/blob/main/results/plots/dendrograms/CitochromeCdendrograms.png)   
#### CytochromeC

Not surprisingly, we can identify as the most informative dendrogram the one calculated using the `ward` merge method, showing more divergent branches and discriminating between two major clusters, consistent with the **silhouette** inference method, but also suggesting 3 as a putative number of clusters. Note that, in contrast to the previous two datasets, this is not the one with the highest agglomerative coefficient (having 0.8447488).

Interestingly, the ones with the highest coefficient (0.9586281) are the least informative trees, those calculated using the `single` merge method.


With these three hierarchical clustering methods, we can compare the differences across different heuristics for resolving the distance information present in our distance matrices. We can note a common tendency in these comparisons: even when different heuristics are used for the agglomerative methods, the topology depends more on the merge method. The `ward` method appears to be more informative in the three sets of proteins, even over the divisive methods. Here, we further break down the information.


## Discussion

Clustering methods are a powerful tool to access analyses of differences in our data, but they are far from being conclusive methods for assessing biological assumptions. They must be used carefully and as complementary tools for confirmation or even hypothesis generation, requiring additional information to reach well-supported conclusions.

In this particular case, the resolution of a dendrogram is not, in any way, a phylogenetic result. Rather, it is a way to explain variation in the distances obtained from normalized scores derived from BLAST alignments. In the cases where the expected clusters were explained (ABC and H1), we can see that choosing between agglomerative heuristics was less relevant than choosing the method by which clusters properly agglomerate. This teaches us that the way in which clusters are considered part of another, or not considered part of the same, is crucial, with methods that account for variation across clusters (such as the `ward` method) being more effective.

Another important point when generating these dendrograms using hierarchical methods is that the agglomerative coefficient only reflects the confidence in the branch structure of the dendrogram, but it does not reflect the confidence we can have in the dendrogram topology itself when discussing biological assumptions.

Finally, in this exercise, the inference method that showed the best performance across all datasets was the `silhouette` method. 
