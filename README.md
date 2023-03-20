# WES-Filtering-Algorithm-for-PKD-Genes-using-ANNOVAR-Outputs

A filtering algorithm for genes of interest to filter down potentially pathogenic variants annotated by ANNOVAR with the aim to find rare potential mutations within genes in ARPKD. 

### Forward:
The script within this GitHub repository was designed and used as part of a PhD thesis to explore rare variants in ARPKD individuals. The script is published here to provide access to the script used within the study for interested individuals. Although the script can filter down to genes of interest, several limitations may limit its use outside this study. However, should any aspects be of interest, it is free to use.     

### Function: 
First, the script extends the available statistics lists by creating new annotations based on the ANNOVAR output. The first new annotation is an average population frequency across the three variant population databases (1,000 Genome Project, NHLBI-ESP project, and Exome Aggregation Consortium). The second annotation is a percentage of the likelihood of being pathogenic based on the predictions across all available prediction programs. Prediction estimates depended on information from dbSNP. The third was a percentage of conserved predictions across all available conservation programs, dependent on the information from dbSNP. The final set of new annotations was pass or fail for the Human Gene Mutation Database (HGMD) and ClinVar databases and CADD scores. All variants were filtered to identify variants located within exon regions or splice sites. A second filter was then applied to isolate rare variants based on the new average population frequency. Any variant with an average population frequency of less than 5% was considered rare for the interests of this analysis. 

The script then splits into modules investigating subsets of genes. First, the module script filters the rare variant subset by gene lists, such as *PKHD1* and *DZIP1L*. After shortlisting the variants by genes of interest, variants were filtered for confirmed pathogenic variants in either ClinVar or HGMD. Variants that failed this filter were then investigated to determine the variant type. Splicing variants, frameshifts, and truncations were then highlighted as potential variants of interest. SNPs were investigated to determine the likelihood of being pathogenic, which was determined as a pathogenic percentage of 30% or higher. Thirty per cent was derived from the need for two or more programs that predicated a likelihood of being pathogenic when four or more predictions were available while allowing variants with less information to also pass through the filter. The final filter detected more novel variants with no pathogenic predictions. Outputs were generated for all rare variants and variants of interest. A final output file of all rare variants that passed the filter parameters and removed the actively investigated genes was also generated as part of the script for further analysis. The script was tested on a small training dataset generated from the sequencing data to ensure that the script could filter and output the data correctly.

![Script Flow Diagram](https://github.com/TpmRichard/WES-Filtering-Algorithm-for-PKD-Genes-using-ANNOVA-Outputs/blob/main/Script%20Flow%20Diagram.png)


### Filter Cutoffs

Currently filter cutoffs are predefined within the code and are outlined below. 

| Application      | Criteria         |
| ---------------- | ---------------- |
| CADD             | Phred Score > 10 |
| FATHMM           | < -1.5           |
| Gerp++           | > 2              |
| LRT              | D                |
| MutationAssessor | H or M           |
| MutationTaster   | A or D           |
| PhyloP           | > 0.5            |
| PolyPhen-2       | > 0.447          |
| SIFT             | < 0.05           |
| SiPhy            | > 2              |

### Limitations:
The script has a few limitations that complicate its application to more general datasets. Many of these aspects relate to a lack of variables within the code. 
1)	Several functions within the current implementation reference specific column names or column positions. Different ANNOVAR outputs may have different names or column positions for similar information. 
2)	The user may wish to use different databases in their studies. This cannot be accounted for in the current implementation.
3)	The user may have different genes of interest that cannot be automated in the current implementation without altering the code. 
4)	Parameters are not as quickly tweaked in the current implementation without altering the code.
5)	Currently, running the code cannot be done by calling a function. 
6)	The script assumes a specific naming convention. 
7)	The script assumes two file inputs: Single Nucleotide Polymorphisms (SNPs) and Insertions and Deletions (INDELs). 

### Acknowledgements

ANNOVAR, which provided the annotations that this script makes use of to highlight variants of interest (Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010). 

Tidyverse package, functions from which this tool makes use. (Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686). 

### Publications:
N/A
