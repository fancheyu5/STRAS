# STRAS
 A snakemake pipeline for genome wild short tandem repeats (STR) annotation and score.
 
1.Set up

First create a conda environment to run this pipeline, by

    conda create -n STRAS
   
    conda activate STRAS
   
    conda install snakemake -c bioconda
   
Then install BEDOPS and r-base, by

    conda install bedops r-base=4 -c bioconda 
   
Install R packages

    conda install r-caret r-tidyverse r-randomForest r-dplyr r-stringr

2.Prepare input .bed

The input file should be a tab-delimited, no header .bed file with 5 columns. The contents of each column are chromosome (without characters ‘chr’), start, end, motif and copy number of this locus. An example is hg002.gangstr.bed . Name your bed with sample1.bed , sample2.bed ...
The reference should be hg38. Please use LiftOver if you use hg19.

3.Start

Edit the first line of snakemake file. Put your sampleIDs in the [] , quot them with "", and seperate by ,

eg. sample=["sample1","sample2"]

Then run

    snakemake

4.Results explain

STRAS could annotate STR loci by genome location, distribution from general population, TAD boundaries, CTCF sites, phenotype of associated gene and a predicted score of being pathogenic. 
STRAS score is blind to loci without population distribution annotation (generated by GangSTR v2.5.0 from 2 1000G cohorts). So GangSTR is recommended for patient STR generation before annotation.             
Pre-mutations could not be classified with high accuracy.
Expansion with peroid length 1 and 2 bp would not be scored. 
Pathogenic insertions could not be identified by STRAS.


5.An example

hg002.gangstr.bed is an example generated from hg002 WGS data. Just run snakemake and it would take a few minutes. Annotated tsv and scored tsv would be in results file.
