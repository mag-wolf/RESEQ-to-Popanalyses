#  RESEQ-to-Popanalyses
-------------------------------

Building a fully automized pipeline for population genomic and conservation genomic analyses.

A main goal is to create high quality SNP files that can be analyzed and visualized in SambaR (https://github.com/mennodejong1986/SambaR).

Can also be used to calculate whole genome statistics like Dxy,HE,pi,tajimasD,wattersonsTheta and test for gene flow using mainly tools
from Simon H Martin (https://github.com/simonhmartin/genomics_general).

Can also be used to model the demographic past of your populations using the MSMC2 framework. This part largely relates to a tutorial
created by Jessica Rick (https://github.com/jessicarick/msmc2_scripts).

Does require perl5, python3 and anaconda environments!

by Magnus Wolf 2023 (magnus.wolf@senckenberg.de)
-------------------------------
M.Sc. Magnus Wolf

PhD-Student

Senckenberg Biodiversity and Climate Research Center (SBiK-F)

Group: Evolutionary Vertebrate Genomics

Georg-Voigt-Straße 14-16, floor 4, room 4.07

60325 Frankfurt am Main, Germany

Installation:
-------------------------------
1.) Download the tarball

2.) copy the tarball and rename it to what you desire

3.) extract the renamed tarball

    tar -xzvf renamed.tar.gz renamed
    
4.) Go into the extracted directory and extract all sub-directories:

    tar -xzvf bin.tar.gz bin
    tar -xzvf RAW_READS.tar.gz RAW_READS

4.) install dependencies via conda:

    conda create --name FASTPenv -c bioconda fastp
    conda create --name MAPPINGenv -c bioconda bwa samtools picard
    conda create --name QUALIMAPenv -c bioconda qualimap
    conda create --name BCFenv -c bioconda bcftools vcftools samtools
    conda create --name PLINKenv -c bioconda plink
    conda create --name PYTHON2env python=2.7
    conda create --name NUMPYenv -c anaconda numpy
    conda create --name MSMC2env -c bioconda msmc2    #This is not a straight forward installation and you might need to manually install dependencies!

Now you are ready to start.

Usage:
-------------------------------
1.) Gather whole genome illumina short read fastq data from all individuals you want to include. Name all input files in a similar
fashion: indivname_1.fq.gz and indivname_2.fq.gz for the forward and reverse reads, respectively. If you have multiple forward
and reverse files, you need to merge these files first. Keep the individual name as short as possible and without special 
characters like " - " or " _ " to avoid that the pipeline messes up names in downstream processes. 

2.) Open the script RESEQ-to-Popanalyses.sh with a text editor of your choice (e.g. nano).

3.) Edit general dependencies. Especially the working directory and the number of threads of parallel processes you wanne use.

4.) Edit the RESEQ-to-Popanalyses scripts and envs if necessary. Usually, if all scripts are at the place where they should be, you shouldn't need to change anything here. 

5.) Edit the popfile.txt and the makeMappabilityMask.py by providing either the names of individuals and their respective populations or the name of the reference fasta-file.

6.) Edit RESEQ-to-Popanalyses parameters and filter settings. Most default setting I used might work for you as well, but some you really need to provide yourself like e.g. maxdp and mindp, the ABBABABA population names and the name of the largest scaffold!

7.) Edit options for RESEQ-to-Popanalyses. Call everything "TRUE" that you want to use. The pipeline 
contains 19 subparts that can be run independently if all other subparts are called "FALSE". By
leaving it as it is, everything will run one by one. However, I wouldn't recommend this as many of these 
steps create large files and will take up more space then you have. Especially domapping, runindivvarcall
and runmsmc2varcall should be handled carefully in this instance!

8.) Now simply run:
    
    bash RESEQ-to-Popanalyses.sh

I suggest piping screen outputs to an error log by adding: 

    2>&1 | tee error.log
 

Here a list of all subparts:
###

rawreadtrimming                   #trimm the provided raw reads and remove adapter using fastp
domapping                         #mapp the trimmed reads to the provided reference using bwa-mem, samtools and picard.
testmapping                       #test the mapping qualities with qualimap
runvariantcalling                 #call variants with your set of mapping files using bcftools mpileup and bcftools call
filtervariances                   #filter the called variances with bcftools filter and vcftools to get the final vcf file as well as a biallelic SNP file
ldpruning                         #prune the SNP file for linkage disequilibrium with bctools prune+
convertvcf                        #convert the SNP file to input files usable for the analysis toolkit SambaR using plink
testkinship                       #test for potential kinships in your set of individuals using plink
thinningforsambar                 #thinn the SNP file further if you have too many
calcdistance                      #calculate raw genetic distances to construct distance based BIONJ trees with SambaR
dodarwindow                       #identify runs of homozygosity to assess inbreeding using DARWINDOW
runindivvarcall                   #call individual variants per mapping file using bcftools mpileup and bcftools call
calchet                           #calculate genome wide heterozygosity using the individual variant files
calcpi                            #calculate the mean pairwise genetic diversity (pi) using the general_genomics tool box from Simon H Martin
calcABBABABA                      #estimate gene flow with a window based ABBA BABA test using the general_genomics tool box from Simon H Martin
angsdtotajima                     #conduct a neutrality test by calculating tajimas D and wattersons theta using ANGSD
dosnpable                         #create a mappability mask using the SNPable tool, used for demographic moddeling with the MSMC2 framework
runmsmc2varcall                   #create individual mask and variance files using bcftools and msmc-tools, used for demographic moddeling with the MSMC2 framework
domsmc2                           #run a demographic moddeling using the MSMC2 framework

###


Good luck!
-------------------------------
