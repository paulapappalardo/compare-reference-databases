# compare-reference-databases

This repo includes R functions to select the best match from BLAST results against different reference databases.

Paula Pappalardo developed the first functions for comparing two specific reference databases (e.g., MIDORI vs. BOLD). Emma Palmer expanded this work, creating general functions to compare up to three reference databases. Emma and Paula work at the Smithsonian Environmental Research Center and have been involved in analyzing metabarcoding data for both the [Coastal Disease Ecology](https://serc.si.edu/labs/coastal-disease-ecology) and the [Marine Invasions laboratories](https://serc.si.edu/labs/marine-invasions-research).

CITATION: Please cite this repo XXX

DEPENDENCIES: The function requires the R packages [dplyr](https://dplyr.tidyverse.org/) and [taxize](https://docs.ropensci.org/taxize/articles/taxize.html).

USAGE: These are R functions taking dataframes or tibbles from blast results as the argument. You need to copy the functions from this repo and load them into your R environment to use them. The specific arguments for each are detailed in the function comments.
