# compare-reference-databases

This repo includes R functions to select the best match from BLAST results against different reference databases.

Paula Pappalardo developed the first functions for comparing two specific reference databases (e.g., MIDORI vs. BOLD). Emma Palmer expanded this work, creating general functions to compare up to three reference databases. Emma and Paula work at the Smithsonian Environmental Research Center, analyzing metabarcoding data for projects in the [Coastal Disease Ecology](https://serc.si.edu/labs/coastal-disease-ecology) and [Marine Invasions](https://serc.si.edu/labs/marine-invasions-research) laboratories.

CITATION: Please cite this repo XXX

DEPENDENCIES: The function requires the R packages [dplyr](https://dplyr.tidyverse.org/), and [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html).

USAGE: These are R functions taking dataframes or tibbles from blast results as the argument. You need to copy the functions from this repo and load them into your R environment to use them. The specific arguments for each are detailed in the function comments.

## Reference databases

For our work, we have been using several reference databases, depending on the project and genetic marker:
* [MIDORI](https://www.reference-midori.info/)
* [BOLD](https://boldsystems.org/data/data-packages/)
* [Lavrador et al. 2023](https://www.mdpi.com/1424-2818/15/2/174)
* [Metazoogene](https://metazoogene.org/mzgdb/)

## Setup

We use this format for BLAST to get results as .tsv files:

``` bash
blastn -task blastn -db dbName -query querySeqsFasta -word_size 11 -max_target_seqs 200 -evalue 0.01 -outfmt "6 qseqid sseqid staxids evalue bitscore length qcovs nident pident" -out blastResults.tsv
```
We use ALL of those columns in our comparing functions, so make sure you have them.

We load results renaming columns to make them more human-readable: 

midori_blast <- fread("blast/blast_results_PWS_to_midori.tsv", fill = T,
                         col.names = c("query_seqid", "result_seqid", "evalue", "bitscore",
                         "length", "percent_coverage", "nident", "percent_identity")) 

We also need to add higher taxonomy for the resulting matches. Our functions require the columns kingdom, phylum, class, order, family, genus, species, and sciname. The column sciname is the taxonomic name of the lowest resolution. For example, if a match was resolved only to the family level, sciname would have the family name, and genus and species would have NA values. 

Finally, we need to keep only the best hit from each database. So for each ASV or OTU, we need only one match per database. For our work, we have been using the best-shared-method from [Pappalardo et al. (2025)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70147).

## Comparison rules

To compare across databases, our functions set up a BLAST percent identity threshold (identity_th) to define what we consider a "good" match. This will depend on the genetic marker and/or taxonomic group so we encourage users to select their own. Comparison across databases were done with the following rules:
* For matches to the same scientific name, we kept the one that had higher similarity measures with the unknown sequence.
* If only one of the matches had a percent identity higher than the BLAST percent identity threshold, we kept that one, regardless of the database of origin.
* When only one database returned a match from BLAST, we kept that one.
* If the scientific names disagree and both matches had a higher than the BLAST percent identity threshold, we kept the one with better taxonomic resolution. If the taxonomic resolution was the same, we would prioritize the local database match, or BOLD over MIDORI.
* If both matches were under the BLAST percent identity threshold, we kept the one with the best similarity measures.


## compare-two-databases.R

This R file includes functions:

* addScinameLevel(): finds the resolution of the scientific name and adds the column _sciname_level_.
* labelFinalMatch(): condenses taxonomy into one string and adds reference database label.
* compareBOLDandMIDORI(): compares matches of BLAST done against BOLD and MIDORI databases.
* pickFinalTax_BoldVsMidori(): keeps only the best match from each database and tags from which database came from.
* compareMLMLandLavrador(): compares matches of BLAST done against MLML and Lavrador databases.
* pickFinalTax_MLMLvsLavrador(): keeps only the best match from each database and tags from which database came from.
* compareGlobalVsCurated(): compares matches from global (midori, bold) to curated (mlml, lavrador) databases.
* pickFinalTax_GlobalVsCurated(): keeps only the best match from each database and tags from which database came from.

Usage:

Here is an example to first compare between two large-scale public databases (bold and midori) and two smaller databases (mlml and lavrador):

```R
bold_vs_midori <- compareBOLDandMIDORI(your_midori_df, your_bold_df, identiy_th) %>%
    pickFinalTax_BoldVsMidori()
    
mlml_vs_lavrador <- compareMLMLandLavrador(your_mlml_df, your_lavrador_df, identiy_th) %>%
    pickFinalTax_MLMLvsLavrador()
```
If you first run only the compare function, you can check the cases where both databases have a very good match (e.g., 99% identity) and whether there are conflicts in the assignments.

You can then compare the outcomes of the two:

```R
your_final_df <- compareGlobalVsCurated() %>%
pickFinalTax_GlobalVsCurated()
```
Now, keep in mind that the taxonomic framework in your databases can be different. It is important to standardize all scientific names to a single taxonomic framework (e.g., GBIF or WoRMS) after you have your final list of matches for each sequence.
