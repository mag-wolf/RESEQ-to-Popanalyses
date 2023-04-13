#!/bin/bash
#
# This is RESEQ-to-Popanalyses, a pipeline created to automatize all my steps for population genomics and conservation genomics.
# In the end you'll have some basic analyses done for your population-genomic data and an input file for SambaR (de Jong et al. 2021).
#
#
# Does require specific conda environments!
#
# by Magnus Wolf 2022 (magnus.wolf@senckenberg.de)
#
# In the following part you need to provide paths and desired settings for one or all of the different scripts you want to run.
#
#################################general dependencies and parameters#####################################################################
source ~/anaconda3/etc/profile.d/conda.sh               #make sure to be able to switch between conda envs
WORKINGDIR=/path/to/working/dir                         #path to the Working directory, without the last #/# symbol!
ulimit -n 2048                                          #get more doable jobs-temp files, I suggest not to change that.
THREADS=N                                               #number of threads
TASKS=n                                                 #number of parallel task to run alignment or trimming processes, keep in mind that during the variant calling step, this part multiply the number of threads used by 10, so 10TASKs with 10THEADS each euals 100 threads used!!!
JAVAMEMORY=NNNG                                         #memory you can afford to provide to java applications like e.g. qualimap. Be careful, java will always take ~10% more than you said it should take. Don't forget the #G# !!!!
REFERENCE=/path/to/reference.fasta                      #path to reference file including the filename, make sure the reference is of good quality,
#                                                       #closely related and not from your population of interest! Also, please name your reference ".fasta" not ".fa" for convinence in some downstream analysis regarding the MSMC2 steps!
#
#
##################################RESEQ-to-Popanalyses scripts and envs#################################################################
RUNMPILEUP=${WORKINGDIR}"/bin/run_mpileup_parallel_v10.sh"            #path to run_mpileup_parallel_v10.sh script for running mpileup variant calling parallel, no need to change if script is in working dir as it should be.
RUNMPILEUPINDIV=${WORKINGDIR}"/bin/run_mpileup_parallel_v10_indiv.sh" #path to run_mpileup_parallel_v10_indiv.sh script that runs variant calling one by one. Only necessary for PSMC analysis!!!!
snpable=${WORKINGDIR}"/bin/seqbility-20091110"                        #path to the snpable script directory. If the dir is in the WorkingDir, no need to change.
CALCDIST=${WORKINGDIR}"/bin/VCF_calcdist.sh"                          #path to the VCF_calcdist.sh script that calculates pairwise raw genetic distances within a multiple-individuals containing vcf file. If the binary is in the WorkingDir, no need to change. 
COUNTHE=${WORKINGDIR}"/bin/countHE.py"                                #path to the countHE.py script to calculate genome wide heterozygosity levels, no need to change if script is in working dir as it should be.
COUNTHEWINDOW=${WORKINGDIR}"/bin/countHEwindow.sh"                    #path to the countHEwindow.py script to calculate heterozygosity per window to identify runs of homozygosity with DARWINDOW, no need to change if script is in working dir as it should be.
DARWINDOWPREP=${WORKINGDIR}"/bin/DARWINDOW_prep.sh"                   #path to the DARWINDOW_prep.sh script to prepare the input-generation of files for the R based DARWINDOW script
DARWINDOWRUN=${WORKINGDIR}"/bin/DARWINDOW_HEperWindow.sh"             #path to the DARWINDOW_HEperWindow.sh script to run the actual input-file generation for the R based DARWINDOW script
BCFTOOLS=~/anaconda3/envs/BCFenv/bin/bcftools                         #path to the bcftools binary, needed for the runvariantcalling and calcdistance step
BCFTOOLSBIN=~/anaconda3/envs/BCFenv/bin                               #path to the bin directory containing bcftools, needed for the dodarwindow step
SAMTOOLS=~/anaconda3/envs/BCFenv/bin/samtools                         #path to the samtools binary, needed for the runvariantcalling step
TABIX=~/anaconda3/envs/BCFenv/bin/tabix                               #path to the samtools binary, needed for the dodarwindow step
VCFTOOLS=~/anaconda3/envs/BCFenv/bin/vcftools                         #path to the samtools binary, needed for the dodarwindow step
BCFTOOLPLUGIN=${WORKINGDIR}"/bin/bcftools_prune/"                     #path to the bcftools plugin directory, if the #bcftools_prune# directory is in the working directory, no need to change this one
SMGENOMICS=${WORKINGDIR}"/bin/genomics_general/"                      #path to the script collection from simon marting calles genomics_general (https://github.com/simonhmartin/genomics_general), no need to change if already in the working directory!
MSMC2TOOLS=${WORKINGDIR}"/bin/msmc-tools/"                            #path to the collection of additional tools provided by the authors of msmc2. If this directory is already in the working dir, no need to change. 
#
FASTPenv=FASTPenv                                       #name of the fastp conda env, should have fastp installed
MAPPINGenv=MAPPINGenv                                   #name of mapping conda env, should have bwa, samtools and picard installed
QUALIMAPenv=QUALIMAPenv                                 #name of qualimap conda env, should have qualimap installed
BCFenv=BCFenv                                           #name of bcftools conda env, should have bcftools, vcftools and samtools installed
PLINKenv=PLINKenv                                       #name of the plink conda env, should have plink installed.
PYTHON2env=PYTHON2env                                   #name of a conda env with python2 installed.
NUMPYenv=NUMPYenv                                       #name of the conda env with numpy installed. Used for every step that involves tools from Simon Martin
MSMC2env=MSMC2env                                       #name of the conda env with msmc2 installed. This is not a straight forward installation and you might need to manually install dependencies like dmd2 and gsl as well as export their binaries to the PATH! Make sure this is running properly!
#
popfile=${WORKINGDIR}"/popfile.txt"                     #For genetic distance, genetic diversity (pi) calculation, ABBA_BABA calculation, and tajimas D estimation, provide a tab seperated file that assignes each of your tested individuals with a respective population name, the individual names should be consistent with your provided short read file names since they will be used throughout the entire script as ID!!!!
makemapmask=${WORKINGDIR}"/makeMappabilityMask.py"      #For the msmc2 analysis, alter this file to specify the path to your MSMC2 directory and your reference name. If you need more information, look at the info text that popps up when you run the dosnpable part for the first time.
#
#
##################################RESEQ-to-Popanalyses parameters and filter settings###################################################
minrawseqlength=40                                      #For raw read filtering: minimal length of reads after trimming, everything shorter will be removed.
qualityphred=15                                         #For raw read filtering: the phred-score quality value that a base is qualified. Default 15 means phred quality >=Q15.
unqualpercent=40                                        #For raw read filtering: how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%
mincomplexity=30                                        #For raw read filtering: threshold for low complexity filter, its default value is 30, which means 30% complexity is required.
minbasequal=13                                          #set filter criteria for variant calling with bcftools mpileup and call (minimal Base quality, default is 13, which roughly corresponds to a probability of 0.05 that base call is incorrect, as given by the formula -10*log10(0.05))
minmapqual=20                                           #set filter criteria for variant calling with bcftools mpileup and call (minimal Mapping quality)
QUICKFILTER=FALSE                                       #if desired, speed up the variant filtering process by filter for every parameter at once. It will trumendely speed up the process but you'll loos the information, which filter setting removed how many sites.
maxdp=1200                                              #set filter criteria for a maximal coverage of SNPs (~3x the total expected coverage with all individual coverages combined)
mindp=100                                               #set filter criteria for a minimal coverage of SNPs (~0.3x the total expected coverage with all individual coverages combined)
minindivcov=3                                           #set a filter criteria for the minimal coverage PER INDIVIDUAL you would like to allow to contribute to variant calling. All others will be masked and thrown out by the maxmissing filtering step. 
maxmissing=0.9                                          #set filter criteria for a maximum of missing data in SNPs with vcftools
LDcutoffr2=0.9                                          #set r2 cutoff for the LD pruning (here a light pruning is used)
LDwindow=1000                                           #set sligins window size for LD pruning
minSNPswanted=1000000                                   #number of desired SNPs for thinning, SambaR running in R studio on a windows pc cant handle more then one million SNPs, hence I used it here as a default. Change if you have a more or less powerfull setup, the SambaR manual on github has more information about performance vs SNPnumber!
ROHwindowsize=10000                                     #For ROH identification with DARWINDOW, set the windows size in bp for which heterozygosity should be calculated, it is usually a good idea to play around with this setting
ROHmincontiglength=3000000                              #For ROH identification with DARWINDOW, set the minimal length of contigs in bp that should be considered for ROH identification, should usually be big enough to find long ROHs of e.g. 1Mbp
minwindowsize=100000                                    #For genetic diversity (pi) calculation, ABBA_BABA calculation and tajimas D estimation that involved a sliding window size, define a window size in bp, using a decent assembly, 100kbp was enough to get ~23k windows to do my calculations on
minsnpsperwindow=100                                    #To avoid using sliding windows without information, setup a threshold for the minimal number of SNPs per window.
ABBABABA_POP1=name1                                     #For ABBA BABA D stats analysis, provide the name of the population that should be placed as P1 in the tested topology: (((P1,P2),P3),O)
ABBABABA_POP2=name2                                     #For ABBA BABA D stats analysis, provide the name of the population that should be placed as P2 in the tested topology: (((P1,P2),P3),O)
ABBABABA_POP3=name3                                     #For ABBA BABA D stats analysis, provide the name of the population that should be placed as P3 in the tested topology: (((P1,P2),P3),O)
ABBABABA_OUT=outname                                    #For ABBA BABA D stats analysis, provide the name of the population that should be placed as O in the tested topology: (((P1,P2),P3),O)
SCAFF1=scaff1                                           #For the MSMC2 analysis, provide the name of your largest scaffold to run the models on. I decided to not run these analyses on the entire genomes due to too high resource requirements and no meaningful differences between complete and 1scaffold predictions. 
#
#
# Now that you provided the dependencies, you can choose which part of which script you wanne run, everything you'll call #TRUE# will run!
# You can find a detailed discription for each step in a makeshift square down below. Additionaly, this discription will pop up if you run the respective
# step.
#
##################################Options for RESEQ-to-Popanalyses#######################################################################
rawreadtrimming=TRUE
domapping=TRUE
testmapping=TRUE
runvariantcalling=TRUE
filtervariances=TRUE
ldpruning=TRUE
convertvcf=TRUE
testkinship=TRUE
thinningforsambar=TRUE
calcdistance=TRUE
dodarwindow=TRUE
runindivvarcall=TRUE
calchet=TRUE
calcpi=TRUE
calcABBABABA=TRUE
angsdtotajima=TRUE
dosnpable=TRUE
runmsmc2varcall=TRUE
domsmc2=TRUE
#
#
#########################################################################################################################################
# Now you're finished. Go into the Working Directory and hit #bash RESEQ-to-Popanalyses.sh#. Then everything you chosed will run one by one.
# I recommend using some sort of log pipeing since long pipelines like this one tend to produce a lot of hidden errors. Try e.g.
# #2>&1 | tee error.log# behind the bash command.
#
# After this line I wouldn't touch anything unless you know what youre doing...
#########################################################################################################################################


if [[ "$rawreadtrimming" = TRUE ]]
        then
        echo "########################start rawreadtrimming######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is rawreadtrimming of RESEQ-to-Popanalyses, a wrapper       #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will trimm all provided fastq.gz files in the    #"
        echo "##RAW_READS# directory using fastp (Shifu Chen et al. 2018). This #"
        echo "#function will basically use all trimming functions provided      #"
        echo "#by fastp and for this to work, make sure to provide cutoffs for  #"
        echo "#minimal read length, phred-score, percentage of unqualified      #"
        echo "#bases and minimal complexity. Default values are provided in the #"
        echo "#raw script. Make sure to only provide paired end data in ONLY    #"
        echo "#two files, merge them beforehand if necessary. Forward files must#"
        echo "#be named with a *1.fq.gz ending while reverse files must be named#"
        echo "#with a *2.fq.gz ending. Will output trimmed files in the         #"
        echo "##TRIMMING# directory.                                            #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir TRIMMING
        cd TRIMMING
        conda activate $FASTPenv
        for file in ./../RAW_READS/*_1.fq.gz;
        do
                bn=`basename $file _1.fq.gz`
                fastp -i "./../RAW_READS/"$bn"_1.fq.gz" -I "./../RAW_READS/"$bn"_2.fq.gz" -o $bn"_trimmed_1.fq.gz" -O $bn"_trimmed_2.fq.gz" -g -3 -l $minrawseqlength -y -Y $mincomplexity -q $qualityphred -u $unqualpercent -c -p -j -h -R -w $THREADS
        done
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$domapping" = TRUE ]]
        then
        echo "########################start domapping############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is domapping of RESEQ-to-Popanalyses, a wrapper             #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will take trimmed paird-end short-read-data and  #"
        echo "#map them onto a provided reference fasta file. The function will #"
        echo "#first attempt to index the provided reference and will output    #"
        echo "#index files at the same position of the reference. Afterwards    #"
        echo "#it will map raw reads with bwa-mem, remove duplicates with       #"
        echo "#picard from the gatk toolkit, replace readgroups with the sample-#"
        echo "#identifier that should be the first part of the fastq-filename,  #"
        echo "#and index resulting files with samtools. Make sure to have all   #"
        echo "#these tools installed via conda in the MAPPINGenv. Will output   #"
        echo "#bam and bai files in the #MAPPING# directory. If everything you  #"
        echo "#wanted went through, I would remove these files as they tend to  #"
        echo "#be pretty large.                                                 #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir MAPPING
        cd MAPPING
        conda activate $MAPPINGenv
        refbasename=`basename $REFERENCE .fasta`
        bwa index -a bwtsw ${REFERENCE}
        samtools faidx ${REFERENCE}
        picard CreateSequenceDictionary R=${REFERENCE} O=$refbasename".dict"
        for file in ./../TRIMMING/*_trimmed_1.fq.gz;
        do
                bn=`basename $file _trimmed_1.fq.gz`
                RGname=$(echo ${bn} | tr "-" "_" | sed s/_//g)
                bwa mem ${REFERENCE} ${WORKINGDIR}"/TRIMMING/"${bn}"_trimmed_1.fq.gz" ${WORKINGDIR}"/TRIMMING/"${bn}"_trimmed_2.fq.gz" -t ${THREADS} | samtools sort -o $RGname".sorted.bam"
                samtools index $RGname".sorted.bam"
                picard MarkDuplicates -Xmx150g MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=$RGname".sorted.bam" O=$RGname".sorted.marked_duplicates.bam" M=marked_dup_metrics.txt REMOVE_DUPLICATES=True
                picard AddOrReplaceReadGroups -Xmx150g I=$RGname".sorted.marked_duplicates.bam" O=$RGname".sorted.rmdup.RG.bam" RGID=$RGname RGPL=illumina RGLB=lib1 RGPU=unit1 RGSM=$RGname
                picard BuildBamIndex -Xmx150g I=$RGname".sorted.rmdup.RG.bam"
                samtools index $RGname".sorted.rmdup.RG.bam"
                rm $RGname".sorted.marked_duplicates.bam"
                rm $RGname".sorted.bam"
                rm $RGname".sorted.bai"
                rm $RGname".sorted.bam.bai"
        done
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$testmapping" = TRUE ]]
        then
        echo "########################start testmapping##########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is testmapping of RESEQ-to-Popanalyses, a wrapper           #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will test the quality of your mapping done in the#"
        echo "#step above via qualimap2 (Okonechnikov et al. 2015). Make sure   #"
        echo "#to have it properly installed in it's own conda env. This step   #"
        echo "#takes rather long and will probably take a lot of memory and     #"
        echo "#threads. It is only worth doing if you're unsure about the       #"
        echo "#quality of your samples or if you are unsure if your reference is#"
        echo "#too unrelated to get good mapping rates. Make also sure to       #"
        echo "#provide how much memory java is allowed to take. Be carefule     #"
        echo "#though, java is a beast and will usually take ~10% more then it  #"
        echo "#should. Will output exhaustive statistics in the #QUALIMAP#      #"
        echo "#directory.                                                       #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir QUALIMAP
        cd QUALIMAP
        conda activate $QUALIMAPenv
        for file in ./../MAPPING/*.sorted.rmdup.RG.bam;
        do
                filename=`basename $file .sorted.rmdup.RG.bam`
                qualimap bamqc -bam $file -outdir "./"$WORKINGDIR"/QUALIMAP/"$filename"_qualimap" -outfile $filename"_qualimap.pdf" -nt $THREADS --java-mem-size=$JAVAMEMORY
        done
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$runvariantcalling" = TRUE ]]
        then
        echo "########################start runvariantcalling####################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is runvariantcalling of RESEQ-to-Popanalyses, a wrapper     #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will run the complete variant calling process by #"
        echo "#utilizing another pipeline written by Menno de Jong, called      #"
        echo "##BAM2VCF_run_mpileup_parallel_HIGHWAY#. For the sake of saving a #"
        echo "#lot of time, the pipeline will attempt to run variant calling in #"
        echo "#parallel on a number of subsets that are defined by the number of#"
        echo "##TASKS# you provided above. Be carful with this setting as       #"
        echo "#you will occupy 10 threads per task and it will get out of hand  #"
        echo "#pretty easily. Variant calling is done by bcftools functions     #"
        echo "#mpileup and call, make sure bcftools and samtools are properly   #"
        echo "#installed in your conda env. The here provided method NEEDS a    #"
        echo "#bcftools version of 1.12 or higher to run, if conda dosen't      #"
        echo "#provide such a recent version, you can specify the path to a     #"
        echo "#manually installed binary above in the user defined seciont. Make#"
        echo "#also sure to provide cutoffs like minimal base-quality and mapp- #"
        echo "#quality. Although, I would probably leave base-quality at the    #"
        echo "#default 13. Eventually, the function will output a vcf file      #"
        echo "#contaning all variances as well as monomorphic sites, ready to be#"
        echo "#filtered in the next step.                                       #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir VARIANCES
        cd VARIANCES
        conda activate $BCFenv
        for file in ./../MAPPING/*.sorted.rmdup.RG.bam;
        do
                ln -s $file
        done
        for file in ./../MAPPING/*.sorted.rmdup.RG.bai;
        do
                ln -s $file
        done
        for file in ./../MAPPING/*.sorted.rmdup.RG.bam.bai;
        do
                ln -s $file
        done
        for file in *.sorted.rmdup.RG.bam;
        do
                echo $file >> bam-list_mpileup.txt
        done
        for file in *.sorted.rmdup.RG.bam;
        do
                bn=`basename $file .sorted.rmdup.RG.bam`
                echo -e ${bn}"\t"${bn} >> popfile.txt
        done
        bash $RUNMPILEUP $TASKS $THREADS $SAMTOOLS $BCFTOOLS $REFERENCE variances bam-list_mpileup.txt popfile.txt $minmapqual $minbasequal
        tabix -p vcf variances.vcf.gz
        wait
        conda deactivate
        rm *my*
        rm *.bcf
        rm *.bam
        rm *.bai
        cd ./../
        fi

if [[ "$filtervariances" = TRUE ]]
        then
        echo "########################start filtervariances######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is filtervariances of RESEQ-to-Popanalyses, a wrapper       #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will take the raw vcf file created above         #"
        echo "#and filter it using the paramaters you provided in the user-     #"
        echo "#defined section. Depending on how sure you are about your        #"
        echo "#used thresholds, you can choose between a QUICKFILTER=TRUE option#"
        echo "#or a QUICKFILTER=FALSE one. The slow one gives you an extra stats#"
        echo "#file to see the impact of your filter settings one by one while  #"
        echo "#the quick option will just run the filter settings alltogether.  #"
        echo "#choosing these parameters can be difficult and usually require   #"
        echo "#gut-feeling. I usually go for a 25% missing data and a 3x as well#"
        echo "#as 0.3x of the expected total coverage as min-depth and max-depth#"
        echo "#which we found after extensive testing to work quite reliable.   #"
        echo "#We also chose against filtering on QUAL (snp-quality) as we found#"
        echo "#that removing uninformative SNPs can bias certain analysis like  #"
        echo "#heterozygosity and SFS calculation.                              #"
        echo "#The results are two new files, a filtered vcf that still contains#"
        echo "#monomorphic sites for e.g. heterozygosity calculation and a      #"
        echo "#filtered vcf file with only biallelic sites for e.g. structure   #"
        echo "#analysis in SAMBAR. Given that you used the slow version,        #"
        echo "#filtering stats are stored in: 'variances_noIndels.vcf.gz'.      #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd VARIANCES
        conda activate $BCFenv
        if [[ "$QUICKFILTER" = TRUE ]]
                then
                bcftools view --threads $THREADS --exclude-types indels variances.vcf.gz -O z > variances_noIndels.vcf.gz
                bcftools filter --threads $THREADS -i "INFO/DP>=$mindp && INFO/DP<=$maxdp" variances_noIndels.vcf.gz -O z > variances_noIndels_DP.vcf.gz
                bcftools filter --threads $THREADS --set-GTs . -e "FMT/DP<"$minindivcov variances_noIndels_DP.vcf.gz -O z -o variances_noIndels_DP_masked.vcf.gz
                vcftools --gzvcf variances_noIndels_DP_masked.vcf.gz --max-missing $maxmissing --recode --recode-INFO-all --stdout | bgzip -c > "variances_noIndels_DP_masked_maxmissing.vcf.gz"
                bcftools view --threads $THREADS --min-alleles 2 --max-alleles 2 variances_noIndels_DP_masked_maxmissing.vcf.gz -O z  > variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz
                rm variances_noIndels.vcf.gz variances_noIndels_DP.vcf.gz variances_noIndels_DP_masked.vcf.gz
        else
                totalvariances=$(zcat variances.vcf.gz | grep -v "#" | wc -l)
                echo "Number of variances in: variances.vcf.gz: "$totalvariances >> variances_filtering-stats.txt
                bcftools view --threads $THREADS --exclude-types indels variances.vcf.gz -O z > variances_noIndels.vcf.gz
                noIndels=$(zcat variances_noIndels.vcf.gz | grep -v "#" | wc -l)
                echo "Number of variances in: variances_noIndels.vcf.gz: "$noIndels >> variances_filtering-stats.txt
                bcftools filter --threads $THREADS -i "INFO/DP>=$mindp && INFO/DP<=$maxdp" variances_noIndels.vcf.gz -O z > variances_noIndels_DP.vcf.gz
                DPfiltered=$(zcat variances_noIndels_DP.vcf.gz | grep -v "#" | wc -l)
                echo "Number of variances in: variances_noIndels_DP.vcf.gz: "$DPfiltered >> variances_filtering-stats.txt
                bcftools filter --threads $THREADS --set-GTs . -e "FMT/DP<"$minindivcov variances_noIndels_DP.vcf.gz -O z -o variances_noIndels_DP_masked.vcf.gz
                vcftools --gzvcf variances_noIndels_DP_masked.vcf.gz --max-missing $maxmissing --recode --recode-INFO-all --stdout | bgzip -c > "variances_noIndels_DP_masked_maxmissing.vcf.gz"
                maxmissinfiltered=$(zcat variances_noIndels_DP_masked_maxmissing.vcf.gz | grep -v "#" | wc -l)
                echo "Number of variances in: variances_noIndels_DP_masked_maxmissing.vcf.gz: "$maxmissinfiltered >> variances_filtering-stats.txt
                bcftools view --threads $THREADS --min-alleles 2 --max-alleles 2 variances_noIndels_DP_masked_maxmissing.vcf.gz -O z  > variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz
                numberofSNPs=$(zcat variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz | grep -v "#" | wc -l)
                echo "Number of variances in: variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz: "$numberofSNPs >> variances_filtering-stats.txt
                rm variances_noIndels.vcf.gz variances_noIndels_DP.vcf.gz variances_noIndels_DP_masked.vcf.gz
                fi
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$ldpruning" = TRUE ]]
        then
        echo "########################start ldpruning############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is ldpruning of RESEQ-to-Popanalyses, a wrapper             #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will take your filtered variance file and will   #"
        echo "#apply an additional filtering for linkage-disequilibrium. LD     #"
        echo "#pruneing has it's ups and downs. One the one hand, you'll remove #"
        echo "#a potential selection bias from similar behaving SNPs in close   #"
        echo "#close proximity, but you'll also artificially remove SNPs and    #"
        echo "#might end up with underepresented regions that still reflect     #"
        echo "#the #TRUE# nature of your genome. Hence, I decided to include    #"
        echo "#this for, but only for, downstream analysis with SAMBAR. The     #"
        echo "#heterozygosity calculation, SFS generation and ROH analysis      #"
        echo "#aren't affected by this step. Also, this step requires a bcftools#"
        echo "#plugin that comes with a manual installation but not in an       #"
        echo "#installation with conda. I included a directory here with the    #"
        echo "#necessary #prune# plugin, but the verion might be different from #"
        echo "#what you installed. Worst case, you need to install bcftools     #"
        echo "#additionally and manually to get a working #plugin# directory. If#"
        echo "#so, don't forget to change the path in the user-defind section.  #"
        echo "For now, this step will only apply a very light pruning due to my #"
        echo "#open doubts about this step...usually juust enough to satisfy    #"
        echo "#reviewers. If you desire a more strict filtering, change the     #"
        echo "#LD r2 cutoff and the window size in the user-defined section.    #"
        echo "#If QUICKFILTER is not TRUE, this step will also include its      #"
        echo "#statistic in the variances_filtering-stats.txt!                  #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd VARIANCES
        conda activate $BCFenv
        export BCFTOOLS_PLUGINS=$BCFTOOLPLUGIN
        bcftools +prune -m $LDcutoffr2 -w $LDwindow variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz -Ov -o variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.vcf.gz
        if [[ "$QUICKFILTER" = TRUE ]]
                then
                echo "skipping stats for ld pruning due to QUICKFILTER=TRUE"
        else
                numberafterLD=$(zcat variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.vcf.gz | grep -v "#" | wc -l)
                echo "Number of variances in: variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.vcf.gz: "$numberafterLD >> variances_filtering-stats.txt
                fi
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$convertvcf" = TRUE ]]
        then
        echo "########################start convertvcf###########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is convertvcf of RESEQ-to-Popanalyses, a wrapper            #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will convert resulting ldpruned vcf files into   #"
        echo "#different files that are used by plink or SAMBAR later on.       #"
        echo "#Requires a working version of PLINK1.9 in the PLINKenv and a     #"
        echo "#working version of vcftools in the BCFenv.                       #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd VARIANCES
        conda activate $BCFenv
        vcftools --gzvcf variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.vcf.gz --plink --out variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned
        conda deactivate
        conda activate $PLINKenv
        plink --file variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned --chr-set 95 --allow-extra-chr --make-bed --recode A --out variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned
        cut -f2 "variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.map" | cut -f1 -d ':' > mycontigs.txt && cut -f2,3,4 "variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.map" > mymap.txt && paste mycontigs.txt mymap.txt > "variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.map" && rm mycontigs.txt mymap.txt
        plink --file variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned --chr-set 95 --allow-extra-chr --make-bed --recode --out variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$testkinship" = TRUE ]]
        then
        echo "########################start testkinship##########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is testkinship of RESEQ-to-Popanalyses, a wrapper           #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will run plink's #--genome# option to create a   #"
        echo "#file that contains Identity-by-descent information. The pi_hat   #"
        echo "#value describes the proportion of sides in IBD, and a reltionship#"
        echo "#of a third degree (cousins) is reached by pi_hat > 0.125, a      #"
        echo "#relationship of a second degree is reached by pi_hat > 0.25. This#"
        echo "#option will not automatically remove these individuals because   #"
        echo "#removing them is pointless for many analyses or easly fixable    #"
        echo "#later on. So this part is optional and more for your interest and#"
        echo "#to explain your results. Requires a working version of PLINK1.9  #"
        echo "#in the PLINKenv!                                                 #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd VARIANCES
        conda activate $PLINKenv
        plink --bfile variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned --genome --chr-set 95 --allow-extra-chr -out variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$thinningforsambar" = TRUE ]]
        then
        echo "########################start thinningforsambar####################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is thinningforsambar of RESEQ-to-Popanalyses, a wrapper     #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will thin your resulting vcf file to end up with #"
        echo "#the number of SNPs you desire. To do so, it will first infere the#"
        echo "#genome size (without repetitive elements) by counting the number #"
        echo "#of sites in the filtered vcf file that still contains monomorphic#"
        echo "#sites. Afterwards it will infere the sliding window size for SNP #"
        echo "#sampling and will repeat a stats step (if QUICKFILTER not TRUE!) #"
        echo "#and a converting thep. This part is necessary, because SambaR    #"
        echo "#needs a lot of ressources when using to many SNPs and more than  #"
        echo "#the default 1Million doesn't really improve SambaRs outcomes.    #"
        echo "#Does require a working BCFenv with bcftools and vcftools         #"
        echo "#installed, as welll as a working PLINKenv with PLINK1.9          #"
        echo "#installed.                                                       #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd VARIANCES
        conda activate $BCFenv
        maxmissinfiltered=$(zcat variances_noIndels_DP_masked_maxmissing.vcf.gz | grep -v "#" | wc -l)
        thinnwindow=$(echo "scale=0; $maxmissinfiltered / $minSNPswanted" | bc)
        vcftools --gzvcf variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned.vcf.gz --thin $thinnwindow --recode --recode-INFO-all --stdout | bgzip -c > "variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned.vcf.gz"
        if [[ "$QUICKFILTER" = TRUE ]]
                then
                echo "skipping stats for thinning due to QUICKFILTER=TRUE"
        else
                numberafterthinn=$(zcat variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned.vcf.gz | grep -v "#" | wc -l)
                echo "Number of variances in: variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned.vcf.gz: "$numberafterthinn >> variances_filtering-stats.txt
                fi
        vcftools --gzvcf variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned.vcf.gz --plink --out variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned
        conda deactivate
        conda activate $PLINKenv
        plink --file variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned --chr-set 95 --allow-extra-chr --make-bed --recode A --out variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned
        cut -f2 "variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned.map" | cut -f1 -d ':' > mycontigs.txt && cut -f2,3,4 "variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned.map" > mymap.txt && paste mycontigs.txt mymap.txt > "variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned.map" && rm mycontigs.txt mymap.txt
        plink --file variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned --chr-set 95 --allow-extra-chr --make-bed --recode --out variances_noIndels_DP_masked_maxmissing_biallelic_LDpruned_thinned
        wait
        conda deactivate
        cd ./../
        fi

if [[ "$calcdistance" = TRUE ]]
        then
        echo "########################start calcdistance#########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is calcdistance of RESEQ-to-Popanalyses, a wrapper          #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will calculate raw genetic distances between     #"
        echo "#every pair possible, given your set of individuals. This         #"
        echo "#calculation is based on a script written by Menno de Jong and    #"
        echo "#requires a correct path defined above, a working BCFenv with     #"
        echo "#bcftools installed and a popfile also provided above. This step  #"
        echo "#is based on the overall filtered variance file that still        #"
        echo "#contains monomorphic sites. Since this file may be gigantic, this#"
        echo "#step might take a loooong time. If your feeling that it takes    #"
        echo "#too long, you may need to thinn the overall (not the SNP) vcf    #"
        echo "#file with the thinningforsambar function and provide this        #"
        echo "#file subsequently to this function here. But this requires you   #"
        echo "#to both times get into the source code and change the respective #"
        echo "#steps yourself. Sorry. Will eventually output a vcfdist.*_*.txt  #"
        echo "#table in the #RAW_DIST# directory. This file can be e.g. used    #"
        echo "#in SambaR to calculate distance based trees and heatmaps!        #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir RAW_DIST
        cd RAW_DIST
        conda activate $BCFenv
        maxindiv=$(cat $popfile | wc -l)
        bash $CALCDIST 1 $maxindiv 1 $maxindiv ${WORKINGDIR}"/VARIANCES/variances_noIndels_DP_masked_maxmissing.vcf.gz" $BCFTOOLS
        conda deactivate
        cd ./../
        fi

if [[ "$dodarwindow" = TRUE ]]
        then
        echo "########################start dodarwindow##########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is dodarwindow of RESEQ-to-Popanalyses, a wrapper           #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will generate input files that can be used in the#"
        echo "#R based runs of homozygosity identification and visualization    #"
        echo "#tool DARWINDOW #https://github.com/mennodejong1986/Darwindow#.   #"
        echo "#Will use the total variance file after filtering and will output #"
        echo "#results into a dedicated directory called #ROH#. Make sure that  #"
        echo "#the paths to the scripts DARWINDOW_prep.sh and                   #"
        echo "#DARWINDOW_HEperWindow.sh are set correctly! make also sure to    #"
        echo "#provide a path to the bcftools bin above. Then you need to       #"
        echo "#define a window size per which HE should be calcultated as well  #"
        echo "#as a minimum size of contigs you want to include in your         #"
        echo "#analysis. Since ROHs can get quite large, I would recommend to   #"
        echo "#use 5Mbp. If your assembly is more fragmented you might wanne    #"
        echo "#go lower. With these files you can now start the actual analysis #"
        echo "#in an R environment, which isn't included in this pipeline.      #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir ROH
        cd ROH
        conda activate $BCFenv
        ln -s ./../VARIANCES/variances_noIndels_DP_masked_maxmissing.vcf.gz
        cp ${WORKINGDIR}"/DARWINDOW_prep.sh" ./
        cp ${WORKINGDIR}"/DARWINDOW_HEperWindow.sh" ./
        bash DARWINDOW_prep.sh $BCFTOOLSBIN variances_noIndels_DP_masked_maxmissing.vcf.gz variances_noIndels_DP_masked_maxmissing 2>&1 | tee DARWINDOW_prep_1.log
        bash DARWINDOW_HEperWindow.sh $BCFTOOLSBIN variances_noIndels_DP_masked_maxmissing.vcf.gz variances_noIndels_DP_masked_maxmissing $ROHwindowsize $ROHmincontiglength 2>&1 | tee DARWINDOW_HEperWindow_1.log
        rm *.sh
        cd ./../
        conda deactivate
        fi

if [[ "$runindivvarcall" = TRUE ]]
        then
        echo "########################start runindivvarcall######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is runindivvarcall of RESEQ-to-Popanalyses, a wrapper       #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will re-run the variant-calling like done above, #"
        echo "#just for each individual individually. This isn't strictly       #"
        echo "#necessary and time and space consuming. However, I noticed that  #"
        echo "#a vcf file called with BCFtools call -m isn't usable for         #"
        echo "#the PSMC pipeline (especially vcf2fq). So, if you want to have   #"
        echo "#a PSMC model, you have to make this extra effort. Be aware though#"
        echo "#that it will take the exact space and time like the intitial     #"
        echo "#variant calling just N x times your samplig size. Will run       #"
        echo "#BCFtools call -c instead of -m. Eventually, all variant files    #"
        echo "#will be directly filtered and the raw variant files will be      #"
        echo "#removed to save space. For that, all filter setting will be      #"
        echo "#re-used like specified for the #filtervariances# section.        #"
        echo "#Be aware that the filtering step will take 5 x your              #"
        echo "#specified number of #TASKS# as threads, so make sure you have    #"
        echo "#enough threads before you start this.                            #"
        echo "#All output files will be stored in the #INDIV_VARIANTS_CALLC#    #"
        echo "#directory.                                                       #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir INDIV_VARIANTS_CALLC
        conda activate $BCFenv
        cd INDIV_VARIANTS_CALLC
        echo "###creating symlinks and tmp files###"
        for file in ./../MAPPING/*.bam;
        do
                ln -s $file
        done
        for file in ./../MAPPING/*.RG.bai;
        do
                ln -s $file
        done
        for file in ./../MAPPING/*.RG.bam.bai;
        do
                ln -s $file
        done
        for file in *.sorted.rmdup.RG.bam;
        do
                bn=`basename $file .sorted.rmdup.RG.bam`
                echo ${bn} >> temp_popfile.txt
        done
        echo "###run variant calling one by one################"
        for file in *.bam;
        do
                filename=$(echo $file | cut -d "." -f1)
                bash $RUNMPILEUPINDIV $TASKS $THREADS $SAMTOOLS $BCFTOOLS $REFERENCE $filename $file $minmapqual $minbasequal
                bcftools index $filename".vcf.gz"
                tabix -p vcf $filename".vcf.gz"
                rm *.bcf
                rm *.my*
        done
        N=$TASKS
        task1(){
                bcftools view --threads 5 --exclude-types indels $line".vcf.gz" -O z > $line"_noIndels.vcf.gz"
        }
        echo "####run indel filtering###"
        for file in temp_popfile.txt;
        do
                while IFS= read -r line;
                do
                        ((i=i%N)); ((i++==0)) && wait
                        task2 "$line" &
                done < "$file"
        done
        wait
        echo "###remove raw variant files###"
        for file in temp_popfile.txt;
        do
                while IFS= read -r line;
                do
                        echo "would have removed: "$line".vcf.gz"
                done < "$file"
        done
        task2(){
                meandepth=$(samtools depth $line".sorted.rmdup.RG.bam" | head -10000000 | awk '{sum += $3} END {print sum / NR}')
                minindivdp=$(echo "scale=4; $meandepth / 4" | bc)
                maxindivdp=$(echo "scale=4; $meandepth * 4" | bc)
                bcftools filter --threads 5 -i "INFO/DP>=$minindivdp && INFO/DP<=$maxindivdp" $line"_noIndels.vcf.gz" -O z > $line"_noIndels_DP.vcf.gz"
        }
        echo "####run coverage filtering, calculates cutoffs by itself by using samtools depth over the first 10Mio sites to save time###"
        for file in temp_popfile.txt;
        do
                while IFS= read -r line;
                do
                        ((i=i%N)); ((i++==0)) && wait
                        task2 "$line" &
                done < "$file"
        done
        wait
        echo " ###remove intermediate indel files###"
        rm *_noIndels.vcf.gz
        task3(){
                vcftools --gzvcf $line"_noIndels_DP.vcf.gz" --max-missing $maxmissing --recode --recode-INFO-all --stdout | bgzip -c > $line"_noIndels_DP_maxmissing.vcf.gz"
        }
        echo "####run maxmissing filtering###"
        for file in temp_popfile.txt;
        do
                while IFS= read -r line;
                do
                        ((i=i%N)); ((i++==0)) && wait
                        task3 "$line" &
                done < "$file"
        done
        wait
        echo "###remove intermediate coverage files###"
        rm *_DP.vcf.gz
        task4(){
                bcftools index $line"_noIndels_DP_maxmissing.vcf.gz"
                tabix -p vcf $line"_noIndels_DP_maxmissing.vcf.gz"
        }
        echo "####run indexing of final vcf files###"
        for file in temp_popfile.txt;
        do
                while IFS= read -r line;
                do
                        ((i=i%N)); ((i++==0)) && wait
                        task4 "$line" &
                done < "$file"
        done
        wait
        echo "###remove tmp files###"
            rm temp_popfile.txt
        rm *.bam
        rm *.bai
        echo "###done###"
        conda deactivate
        cd ./../
        fi

if [[ "$calchet" = TRUE ]]
        then
        echo "########################start calchet##############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is calchet of RESEQ-to-Popanalyses, a wrapper               #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will calculate genome wide heterozygosity for    #"
        echo "#all individuals using the indiviual variant files created in the #"
        echo "#step above. It will go through each file and will count sites    #"
        echo "#with unequal reference and alternative alleles. In the end it    #"
        echo "#will calculate heterozygosity by divinding the number of         #"
        echo "#potential heterozygous sites through the total number of sites.  #"
        echo "#Will output resulting statistics in a #het_stats_indiv.txt# file.#"
        echo "#Make sure that the path in the user-defind section for the script#"
        echo "#countHE.py is set correctly and in the WorkingDir.               #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir DIVER_STATS
        cd INDIV_VARIANTS_CALLC
        task5(){
                python $COUNTHE $file >> ./../DIVER_STATS/HE_stats_indiv.txt
        }
        N=$TASKS
        for file in *maxmissing.vcf.gz;
        do
                ((i=i%N)); ((i++==0)) && wait
                task5 "$file" &
        done
        wait
        cd ./../
        fi

if [[ "$calcpi" = TRUE ]]
        then
        echo "########################start calcpi###############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is calcpi of RESEQ-to-Popanalyses, a wrapper                #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will calculate the mean pairwise genetic (pi)    #"
        echo "#diversity for individuals of every populations in your provided  #"
        echo "#popfile.txt. Make sure that this file containts a tab seperated  #"
        echo "#list of all your individuals with an assigned population name    #"
        echo "#in the second column. Individual names should be consistend with #"
        echo "#the raw seuquence file names and hence read groups. Will use     #"
        echo "#two script from the population genetic tool collection hosted by #"
        echo "#simon martins (GitHub:simonhmartin/genomics_general), namely     #"
        echo "#parseVCF.py to convert the main vcf file to a #geno.gz# file and #"
        echo "#popgenWindows.py for the actual window based assesment of pi.    #"
        echo "#make sure to have a working conda env with numpy installed as    #"
        echo "#well as a correct path to the #/genomics_general# directory.     #"
        echo "#Eventually, results can be found in the #DIVER_STATS# directory  #"
        echo "#in a file denoted as: diver_stats_window_*.csv.gz.               #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir DIVER_STATS
        cd DIVER_STATS
        conda activate $NUMPYenv
        cat $popfile | cut -f2 | sort | uniq > temp_pops.txt
        for file in temp_pops.txt;
        do
                while IFS= read -r line;
                do
                        echo "_" >> temp2_pops.txt
                        echo $line >> temp2_pops.txt
                done < "$file"
        done
        cat temp2_pops.txt | tr "\n" "_" | sed "s/__/ -p /g" | sed "s/_//g"  > temp3_pops.txt
        tempcmd=$(cat temp3_pops.txt)
        rm temp_pops.txt temp2_pops.txt temp3_pops.txt
        ln -s ${WORKINGDIR}"/VARIANCES/variances_noIndels_DP_masked_maxmissing.vcf.gz"
        python $SMGENOMICS"/VCF_processing/parseVCF.py" -i "variances_noIndels_DP_masked_maxmissing.vcf.gz" --skipIndels | bgzip > "variances_noIndels_DP_masked_maxmissing.geno.gz"
        python $SMGENOMICS"/popgenWindows.py" -w $minwindowsize -m $minsnpsperwindow -s $minwindowsize -g "variances_noIndels_DP_masked_maxmissing.geno.gz" -o "diver_stats_window_"$minwindowsize".csv.gz" -f phased -T $THREADS""$tempcmd --popsFile $popfile --minData 0.5 --writeFailedWindows
        conda deactivate
        cd ..
        fi

if [[ "$calcABBABABA" = TRUE ]]
        then
        echo "########################start calcABBABABA#########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is calcABBABABA of RESEQ-to-Popanalyses, a wrapper          #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will run a simple ABBA BABA test using tools from#"
        echo "#the population genetics tool collection hosted by simon martins  #"
        echo "#(GitHub:simonhmartin/genomics_general), namely #parseVCF.py# to  #"
        echo "#convert the main SNP file to a #geno.gz# file and                #"
        echo "##ABBABABAwindows.py# to conduct a sliding-window based D-stat    #"
        echo "#gene flow test. Make sure to provide a tab seperated population- #"
        echo "#file that denotes the individual name and the respective         #"
        echo "#population assignments. Also, you need to define in the user-    #"
        echo "#section, which population should take which place in the from    #"
        echo "#the tool tested topology: (((P1,P2),P3),O). If unsure about the  #"
        echo "#true topology, I would recommend to test all combinations of     #"
        echo "#interest and report all results. Make sure to have a working     #"
        echo "#conda env with numpy installed, as well as a correct file to the #"
        echo "##/genomics_general# directory in the user section. Will output   #"
        echo "#results to a #ABBABABA_*.output.csv# table in a dedicated dir    #"
        echo "#called #ABBA_BABA#.                                              #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir ABBA_BABA
        cd ABBA_BABA
        conda activate $NUMPYenv
        ln -s ${WORKINGDIR}"/VARIANCES/variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz"
        python $SMGENOMICS"/VCF_processing/parseVCF.py" -i "variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz" --skipIndels | bgzip > "variances_noIndels_DP_masked_maxmissing_biallelic.geno.gz"
        python $SMGENOMICS"/ABBABABAwindows.py" -g "variances_noIndels_DP_masked_maxmissing_biallelic.geno.gz" -f phased -o "ABBABABA_P1"$ABBABABA_POP1"_P2"$ABBABABA_POP2"_P3"$ABBABABA_POP3"_O"$ABBABABA_OUT".output.csv" -w $minwindowsize -m $minsnpsperwindow -s $minwindowsize -P1 $ABBABABA_POP1 -P2 $ABBABABA_POP2 -P3 $ABBABABA_POP3 -O $ABBABABA_OUT -T $THREADS --minData 0.5 --popsFile $popfile --writeFailedWindows
        rm "variances_noIndels_DP_masked_maxmissing_biallelic.vcf.gz"
        conda deactivate
        cd ..
        fi

if [[ "$angsdtotajima" = TRUE ]]
        then
        echo "########################start angsdtotajima#########################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is angsdtotajima of RESEQ-to-Popanalyses, a wrapper         #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will conduct an estimation of tajimas D, a       #"
        echo "#measure that tells you if a population evolves neutrally or not. #"
        echo "##very simplyfied, inform yourself!!!!# Uses ANGSD and generates  #"
        echo "#a folded site frequency spectrum first that is then used to      #"
        echo "#calculate theta and tajimas D per populations. Make sure         #"
        echo "#to have a working conda env with ANGSD installed. Especially     #"
        echo "#the first step that calculates site allele frequency likelihoods #"
        echo "##saf# takes an extreme amount of RAM! For my largest population  #"
        echo "#of 14 individuals I used some 300-400Gb! Make sure your system   #"
        echo "#supports these requirements. Eventually, this step will conduct  #"
        echo "#an overall and a sliding window based estimation of wattersons   #"
        echo "#theta, tajimas D and a courple of other useful statistics that   #"
        echo "#can be used to discuss e.g. selection and recent bottlenecks.    #"
        echo "#Will run this analysis on every mentioned population in your     #"
        echo "#population file provided above. Also, make sure that the         #"
        echo "#individual names in the population file and the basenames of the #"
        echo "#mapping files in the MAPPING directory are identical! For        #"
        echo "#more information how you get from a side frequency spectrum to   #"
        echo "#tajimas D, see http://popgen.dk/angsd/index.php/Thetas,Tajima    #"
        echo "#,Neutrality_tests and the papers cited there!                    #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir NEUTRAL_STATS
        cd NEUTRAL_STATS
        conda activate $ANGSDenv
        cat $popfile | cut -f2 | sort | uniq > temp_pops.txt
        for file in temp_pops.txt;
        do
                while IFS= read -r line;
                do
                        for file2 in $popfile;
                        do
                                while IFS= read -r line2; #just to create a bam-list per population
                                do
                                        indiv=$(echo $line2 | cut -f1)
                                        pops=$(echo $line2 | cut -f2)
                                        if [[ "$pops" == "$line" ]]
                                        then
                                                echo ${WORKINGDIR}}"/MAPPING/"$indiv".sorted.rmdup.RG.bam" >> "temp_"$line"_bam-list.txt"
                                        fi
                                done < "$file2"
                        done
                        echo "run angsd doSaf on "$line
                        angsd -bam "temp_"$line"_bam-list.txt" -doSaf 1 -anc $REFERENCE -GL 1 -nThreads $THREADS -P $THREADS -out $line"SFS"
                        echo "run realSFS on "$line
                        realSFS $line"SFS.saf.idx" -P $THREADS -fold 1 > $line".fsfs"
                        echo "run saf2theta on "$line
                        realSFS saf2theta $line"SFS.saf.idx" -sfs $line".fsfs" -outname $line"_ftheta" -P $THREADS
                        echo "run thetaStat do_stat total on "$line
                        thetaStat do_stat $line"_ftheta.thetas.idx" -outnames $line"_thetasingle.theta.gz"
                        echo "run thetaStat do_stat window-based on "$line
                        thetaStat do_stat $line"_ftheta.thetas.idx" -win $minwindowsize -step $minwindowsize -outnames $line"_thetawindow.theta.thetasWindow.gz"
                done < "$file"
        done
        echo "finished angsdtajima run, removing temp files"
        rm temp_pops.txt temp*_bam-list.txt
        cd ..
        fi

if [[ "$dosnpable" = TRUE ]]
        then
        echo "########################start dosnpable############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is dosnpable of RESEQ-to-Popanalyses, a wrapper             #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will create a mappability mask, necessary to     #"
        echo "#run a MSMC2 analysis. It will take your specified reference      #"
        echo "#assembly.fasta file and cut it in 35-mers and mapp them onto     #"
        echo "#the reference file again using bwa. Make sure you have a working #"
        echo "#mapping conda env like described above. Afterwards, some format  #"
        echo "#steps are done using Heng Li's scripts. Make sure that the path  #"
        echo "#to the #seqbility-20091110# tool collection is correct. For more #"
        echo "#infos, see #https://lh3lh3.users.sourceforge.net/snpable.shtml#. #"
        echo "#Eventually, the mask is created using a pyton script from the    #"
        echo "#authors of MSMC2 that you need to alter first!!!! To do this,    #"
        echo "#type #nano makeMappabilityMask.py# and change line 26 and 30.    #"
        echo "#More specifically, in line 26, you need to change the            #"
        echo "##/path/to/the/WorkingDir# to the path of your working dir and    #"
        echo "#replace #referencename# with the name of your reference assembly #"
        echo "#fasta file without the #.fasta# ending. In line 30 its pretty    #"
        echo "#similar. Replace the path to your working dir and the reference  #"
        echo "#name.                                                            #"
        echo "#Apart of this, this step is incredible space costly and will     #"
        echo "#create an enormous amounts of tempory files that will be         #"
        echo "#removed to an extend (will run #rm x*#) but while running this,  #"
        echo "#my working dir got dumped with some extra ~600Gb using a         #"
        echo "#2.4Gbp large reference genome. Make sure you have this space!    #"
        echo "#Make also sure to have a conda env with Python2 running!         #"
        echo "#                                                                 #"
        echo "#Eventually, you'll have a #referencename_masked_mask.35.50.fa    #"
        echo "#file in your #/MSMC2/MASK/# directory ready for the next step.   #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        mkdir MSMC2
        cd MSMC2
        mkdir MASK
        cd MASK
        bn=`basename $REFERENCE .fasta`
        echo "start creating 35-mers from reference fasta..."
        $snpable"/splitfa" $REFERENCE 35 | split -l 20000000
        conda activate $MAPPINGenv
        echo "start mapping 35-mers to reference fasta..."
        bwa aln -t $THREADS -R 1000000 -O 3 -E 3 $REFERENCE $bn"_split.35" > $bn"_split.35.sai"
        bwa samse -f $bn"_split.35.sam" $REFERENCE  $bn"_split.35.sai" $bn"_split.35"
        conda deactivate
        echo "start generating snpable mask..."
        $snpable"/gen_raw_mask.pl" $bn"_split.35.sam" > $bn"_rawMask.35.fa"
        $snpable"/gen_mask" -l 35 -r 0.5 $bn"_rawMask.35.fa" > $bn"_mask.35.50.fa"
        conda activate $PYTHON2env
        echo "start creating msmc2 mask..., dont forget to change the makeMappabilityMask.py script!!!"
        python $makemapmask
        conda deactivate
        rm x* *.sam *.35
        cd ..
        cd ..
        fi

if [[ "$runmsmc2varcall" = TRUE ]]
        then
        echo "########################start runmsmc2varcall######################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is runmsmc2varcall of RESEQ-to-Popanalyses, a wrapper       #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will create the vcf files and individual masking #"
        echo "#files necessary to run a proper msmc2 estimation. Make sure to   #"
        echo "#have a running BCFenv with bcftools and samtools installed. Will #"
        echo "#also require the tool collection #msmc-tools# from the authors   #"
        echo "#of msmc2, usually provided with this pipeline. Make sure the     #"
        echo "#paths are correct. You also need to provide the name of your     #"
        echo "#preferably largest scaffold in the #SCAFF1# variable above. After#"
        echo "#a lot of testing, I decided to not run this analysis on the      #"
        echo "#entire genomes because of the massive resource requirements      #"
        echo "#without visible improvments. So, MSMC2 is only run on the largest#"
        echo "#scaffold. Will first go through your mapping files in the        #"
        echo "##MAPPING# dir and calculate the mean depth per mapping file.     #"
        echo "#Then will run bcftools mpileup and call -c and pipe the output   #"
        echo "#together with the mean depth to #bamCaller.py# from msmc-tools.  #"
        echo "#Also, mapping and basequal cutoffs are used like in the          #"
        echo "#above variant calling steps.                                     #"
        echo "#Will output both the vcf.gz and the bed.gz files into the        #"
        echo "##MSMC2/VCFandBED# directory. These files are eventually used in  #"
        echo "#the next step, the msmc2 run.                                    #"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd MSMC2
        mkdir VCFandBED
        cd VCFandBED
        conda activate $BCFenv
        task6(){
                filename=$(echo $file | cut -d "." -f1)
                meandepth=$(samtools depth -r $SCAFF1 $file | head -3000000 | awk '{sum += $3} END {print sum / NR}' | cut -d "." -f1)
                bcftools mpileup -A --min-MQ $minmapqual --min-BQ $minbasequal -C50 -a "DP,AD" -r $SCAFF1 -r $REFERENCE $file | bcftools call -c -V indels | $MSMC2TOOLS"/bamCaller.py" $meandepth $filename"."$SCAFF1".bed.gz" | gzip -c > $filename"."$SCAFF1".vcf.gz"
        }
        N=$TASKS
        for file in ${WORKINGDIR}"/MAPPING/*.sorted.rmdup.RG.bam;
        do
                ((i=i%N)); ((i++==0)) && wait
                task6 "$file" &
        done
        wait
        cd ..
        cd ..
        fi

if [[ "$domsmc2" = TRUE ]]
        then
        echo "########################start domsmc2##############################"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "#This is domsmc2 of RESEQ-to-Popanalyses, a wrapper               #"
        echo "#function for all kinds of population genetic and conservation    #"
        echo "#genetic analyses that can be done using paired-end re-sequencing #"
        echo "#data, this part will create final input files using the vcf      #"
        echo "#and bed files created above and will subsequently run the actual #"
        echo "#msmc2 analysis. Before you start: after a lot of testing I       #"
        echo "#figured that running msmc2 for entire populations did only result#"
        echo "#in completely equal demographic trajectories for all my whale    #"
        echo "#species. To find some differences, I needed to run this analysis #"
        echo "#per individual which makes this here a close resemblence of a    #"
        echo "#PSMC anlysis. Also, I noticed that running this thing with       #"
        echo "#entire genomes makes this step incredible resource intensive     #"
        echo "#and I needed some 500GB of RAM for 5ish parallel individual runs,#"
        echo "#without it improving the graphs compared to runs on only one     #"
        echo "#large scaffold. Also, the installation of MSMC2 via conda is     #"
        echo "#not straight forward (although this might be fixed while you're  #"
        echo "#reading this) and I needed to manually install gsl and dmd2      #"
        echo "#as well as exporting their bins to the PATH variable. Make sure  #"
        echo "#that msmc2 is running properly after installing it with conda!   #"
        echo "#The first step also requires a correct path to the msmc-tool     #"
        echo "#collection. Make sure that there is only one #*.mask.bed.gz# file#"
        echo "#in your #MSMC2/MASK# directory. I have no idea how you gonna     #"
        echo "#call this file so it relies on the occurence of only one file    #"
        echo "#with this specific ending. In case of running msmc2 on only one  #"
        echo "#chromosome, the authors of msmc2 suggest to define the           #"
        echo "#coalescent units as follows: -p 1*2+15*1+1*2. Which is now the   #"
        echo "#default mode of this script. To change this, you need to go      #"
        echo "#into the source code of this script, sorry. Will also run 100    #"
        echo "#baumwelch iterations and no bootstrap replications since I had   #"
        echo "#actual replicates. If you want to play with the code,            #"
        echo "#bootstraps and population-wide estimations, I would recommend    #"
        echo "#having a look at #GitHUB:jessicarick/msmc2_scripts# tutorial!    #"
        echo "#You can find the output tables that you need for plotting        #"
        echo "#with e.g. ggplots geom_stairs function in the #MSMC2_OUTPUT# dir.#"
        echo "#                                                                 #"
        echo "#Good luck!                                                       #"
        echo "#                                                                 #"
        echo "#                                                                 #"
        echo "###################################################################"
        cd MSMC2
        mkdir MULTIHETSEP
        cd MULTIHETSEP
        task7(){
                bn=`basename $file .vcf.gz`
                indivname=$(echo $bn | cut -d "." -f1)
                $MSMC2TOOLS"/generate_multihetsep.py" --mask=${WORKINGDIR}"/MSMC2/VCFandBED/"$indivname"."$SCAFF1".bed.gz" --mask=${WORKINGDIR}"/MSMC2/MASK/"*".mask.bed.gz" $file > ${WORKINGDIR}"/MSMC2/MULTIHETSEP/"$indivname"."$SCAFF1".multihetsep.txt"
        }
        N=$TASKS
        for file in ${WORKINGDIR}"/MSMC2/VCFandBED/"*".vcf.gz";
        do
                ((i=i%N)); ((i++==0)) && wait
                task7 "$file" &
        done
        wait
        cd ..
        mkdir MSMC2_OUTPUT
        conda activate $MSMC2env
        task8(){
                bn=`basename $file2 .multihetsep.txt`
                msmc2 -t $THREADS -p 1*2+15*1+1*2 -i 100 -o ${WORKINGDIR}"/MSMC2/MSMC2_OUTPUT/"$bn".msmc2" -I 0,1 $file2
        }
        for file2 in ${WORKINGDIR}"/MSMC2/MULTIHETSEP/"*".multihetsep.txt";
        do
                ((i=i%N)); ((i++==0)) && wait
                task8 "$file2" &
        done
        wait
        conda deactivate
        cd ..
        fi

##########################################################################################################################################


#What are you doing down here? The script is over bud!
