# Practical 1 - part 2: Repeat annotation using uncurated and curated libraries

## Introduction

In this part of the tutorial we learn how to run `RepeatMasker` using uncurated and curated custom libraries. 

We're going to mask the same tardigrade genome with two different TE libraries in order to compare how curation can change (and improve) annotation. We will learn how to curate a library in `Practical 2` but for now, just look at the results of such process.

The library and results of this tutorial are taken from the efforts done in previous editions of this course and now published in [Teaching transposon classification as a means to crowd source the curation of repeat annotation – a tardigrade perspective](https://link.springer.com/article/10.1186/s13100-024-00319-8) Peona et al. 2024.

We mask:

- `ramVar` genome with `ramVar_rm1.0.fasta` that we generated in `Practical 1 -part 1`
- `ramVar` genome with `ramVar_rm2.0.fasta` which is a (partially) curated transposable element library

After masking the genomes, we are going to visualise the differences in annotation by creating repeat landscapes.

## RepeatMasker

Let's give a look at the different options by typing:

```bash
conda activate te25
RepeatMasker
#or for having a more detailed help
RepeatMasker -help
```

Since the manual is quite long, you can temporarly save it in a file

```bash
RepeatMasker -help > RMSK_man
less RMSK_man
```

Some important options you should look at are

```
 # number of cores to use
 -pa(rallel) [number]

 # options for the speed of the repeat annotation
 -s 	slow *RECOMMENDED*
 -q 	quick
 -qq 	very quick

 # options for including or excluding some repeat categories
 -nolow 	do not mask low complexity
 -noint 	do not mask interspersed repeats - mask only simple/low complexity
 -norna 	do not mask RNA genes
 -alu 		mask only Alu (relevant for primates)
 -div 		mask repeats below a certain divergence threshold (e.g., useful when you want to know only about young repeats)

# option for using a specific repeat library!!
-lib 	specify the path to your custom library *RECOMMENDED*

# option to specify a species name
-species 	RepeatMasker will use consensus sequences specific to that species (e.g., if specify mouse, it will use sequences specific to mouse, sequences present in all Rodents and in all Mammals). -lib and -species are mutually exclusive *RECOMMENDED*

# option for dealing with contamination (bacterial insertion sequences)
-is_only
-is_clip
-no_is

# options for masking
-gccalc 	calculates the GC content of the sequence to better distinguish repeats from low complexity sequences *RECOMMENDED*

# options for the output
-dir 		you can specify an output folder *RECOMMENDED*
-a 		write alignment file (.align) *RECOMMENDED*
-xsmall		returns a fasta file with soft-masked repeats *RECOMMENDED*
-x 		returns a fasta file with hard-masked repeats
-html 		returns the RMSK output as an additional html file
-gff 		returns the RMSK output as an additional GFF file *IMPORTANT*
-e(xcln) 	outputs a table file with a summary of the repeat proportions *RECOMMENDED*
```

:question: Which options would you use?

Here I show you how to run RepeatMasker BUT, as before, we do not have the time to wait for the software to finish.

Start masking the `ramVar` genome:

```bash
# create working directories
mkdir -p ~/TE25/Practical1/RMSK/ramVar_rm1.0  ~/TE25/Practical1/RMSK/ramVar_rm2.0

# copy the required libraries
cp ~/Share/TE25/Practical1/Data/ramVar_rm* ~/TE25/Practical1/Data/

# be sure to be in the right directory
cd ~/TE25/Practical1

# run RepeatMasker with uncurated ramVar_rm1.0.fasta library
RepeatMasker -pa 1 -a -xsmall -gccalc -excln -gff -s -dir RMSK/ramVar_rm1.0 -lib Data/ramVar_rm1.0.fasta Data/ramVar.fasta
```

You can now run `RepeatMasker` on `ramVar` genome using the (partially) curated library `ramVar_rm2.0.fasta`.

I put all the output files in `~/Share/TE25/Practical1/RMSK`. You can copy everything in your working directory.

```bash
cp -r ~/TE25/Share/Practical1/RMSK ~/TE25/Practical1/RMSK
```

In the output directories then you will find several output files:

#### .out
the main and most important output file of RepeatMasker. This files contains 15 columns describing in detail each repeat found in the genome.

Which is the column containing the divergence from consensus?
What do the columns 8 and 14 mean?

#### .align
the `.align` file contains all the alignments found by RepeatMasker and give some supplementary information about the alignment statistics. This file is important to keep for some downstream analysis below.

#### .tbl
this is a summary table (the one produced thanks to the option `-excln`) containg useful information about the portion of sequence masked and proportion of repeats.

#### .masked
this is the soft-masked version of the input fasta file

#### .cat.gz
this is the version of the `.align` file before the `ProcessRepeats` step of RepeatMasker

Now we explore some useful tools and commands to investigate and manipulate the RepeatMasker output files. Many of the scripts listed below are found in the `util` folder within the installation folder of RepeatMasker itself.

## Create landscape plots

### Re-calculate the divergence adjusted for GC content and create the TE landscapes

The presence of methylated cytosines at CpG sites (in animals) often causes the hypermutability of Cs into Ts (leading to TpG and CpA sites) and this mechanism can cause an overestimation of the percentage of divergence from the consensus sequence. Therefore, there is a Perl scripts in the utils of RepeatMasker that can re-calculate the divergence percentages taking into account the presence of CpG sites starting from the `.align` file.

**Note** the CpG re-calculation given here is heavily biased for animals and other organisms in which the context CpG is methylated so if you are using your own data and you know CpGs are not methylated you can skip this step. Also it should be possible to modify the script to consider other methylation contexts.

Here 

```bash
cd ~/TE25/Practical1/RMSK/ramVar_rm1.0

# Calculate the divergence with Kimura 2-parameter distances (excluding CpG sites)
calcDivergenceFromAlign.pl -s ramVar.fasta.align.divsum -a ramVar.fasta.align_with_div ramVar.fasta.align

cd ../ramVar_rm2.0

# Calculate the divergence with Kimura 2-parameter distances (excluding CpG sites)
calcDivergenceFromAlign.pl -s ramVar.fasta.align.divsum -a ramVar.fasta.align_with_div ramVar.fasta.align
```

What is `-s`? What is `-a`?

Some more info about the CpG sites directly from the help of the `calcDivergenceFromAlign.pl` script:

>Treat `CG` dinucleotide sites in the consensus sequence as follows:
>Two transition mutations are counted as a single transition, one transition is counted as 1/10 of a standard transition, and transversions are counted normally (as the would outside of a CpG site). This modification to the Kimura 2 parameter model accounts for the extremely high rate of mutations in at a CpG locus.

#### :warning: Note
The `.divsum` file is very useful to prioritise the consensus sequences to manually curate. You can sort the library for number of bases annotated by each consensus sequence, or sort it by mean divergence so you can decide to first curate the most abundant and/or the youngest ones.

### Make a repeat landscape plot

To make a nice and fast landscape plot we can easily use another RepeatMasker script that will take the output from the previous step as input.

```bash
# Make repeat landscape from CpG-corrected RepeatMasker .align file and average divergence file
# it requires the size of the genome

cd ~/TE25/Practical1/RMSK/ramVar_rm1.0
SIZE=`grep -i 'total length' ramVar.fasta.tbl | tr -s ' ' | cut -f 3 -d ' '`
createRepeatLandscape.pl -div ramVar.fasta.align.divsum -g $SIZE > ramVar.fasta.align.divsum.html

cd ../ramVar_rm2.0
SIZE=`grep -i 'total length' ramVar.fasta.tbl | tr -s ' ' | cut -f 3 -d ' '`
createRepeatLandscape.pl -div ramVar.fasta.align.divsum -g $SIZE > ramVar.fasta.align.divsum.html
```

What is `-g`?
The number given in the command line is specific for the ramVar genome, when running the command on another genome you must change that parameter according to the genome size! You can find the size of the fasta analysed in the `.tbl` file

Download locally the `.html` file and give it a look!

I prepared all the `.html` files for you, compare them all!

Do they look the same? What did it change?

For more comparisons and more details, you can look at this comprehensive [figure](https://link.springer.com/article/10.1186/s13100-024-00319-8#Fig4): *Ramazzottius varieornatus* is on the bottom row.

![landscape_plot](figure4.png)

--- 

## OPTIONAL

### Convert an align file into .out

From the `calcDivergenceFromAlign.pl`, we obtained a new `.align` file with the CpG-corrected divergence but it would be nice to have it in the form of the `.out` file again so here's a one-liner awk command to do so:

```bash
awk '{if(index($12, "(") != 0){print $1,$2,$3,$4,$5,$6,$7,$7-($6-1),$8,"+",$9,$10,$11,$12,$13,$14}else if($9 ~ /C/){print $1,$2,$3,$4,$5,$6,$7,$7-($6-1),$8,$9,$10,$11,$12,$13,$14,$15}}' OFS="\t" ramVar.fasta.align_with_div | sed 's/#/\t/' > ramVar.fasta.align_with_div.noCpG.size.out
```

### Convert the .out file into a BED file

Sometimes it's very handy to have a BED file to use with BEDTools (e.g., intersect repeat annotation with other types of annotations; Quilan and Hall 2010).
We show you three ways to convert the `.out` into BED

1. Make 0-formatted BED file from 1-formatted RepeatMasker .out file (sequence, masking start, masking end) with awk:

```bash
awk '(NR>3){$6=$6-1; print $5,$6,$7}' OFS="\t" ramVar.fasta.out > ramVar.fasta.out.bed
```

2. Make BED6 file from RepeatMasker `.out` file (BED file with TE information):

```bash
awk '(NR>3){$6=$6-1; print $5,$6,$7,$10,$11,$9}' OFS="\t" ramVar.fasta.out > ramVar.fasta.out.bed6
```

### Convert the .out file into a GFF file

GFF files are as useful as BED files therefore let's see how to convert the `.out` file into a GFF using a RepeatMasker util script:

```bash
rmOutToGFF3.pl  ramVar.fasta.out > ramVar.fasta.out.gff
```

Note that you can always tell RepeatMasker itself to generate a GFF file in addition to the out file with the option `-gff`.

### Get sequences from the .out file 

The script `rmOut2Fasta.pl` is an easy and fast way to get a fasta file from the `.out` file:

```bash
rmOut2Fasta.pl -fasta ramVar.fasta -out ramVar.fasta.out
```

It can also be done using a BED file (+ the original fasta file) and the function `bedtools getfasta` from BEDTools.

### Make a hard-masked (NNN) version of the soft-masked (lowercase) .masked file 

Sometimes it is useful to have a hard masked version of the genome and we can do this with a simple perl command:

```bash
perl -e 'while(<>) { if ($_ =~ /^>.*/) { print $_; } else { $_ =~ tr/acgt/N/; print $_;}}' < ramVar.fasta.masked > ramVar.fasta.masked.hard
```

You can also tune the soft and hard masking of a genome by using a custom BED track of repeats and `bedtools maskfasta`.

### Programs in RepeatMasker for downstream analyses:
There are many other scripts included in the RepeatMasker installation that are useful for downstream analysis.

**Scripts in the main RepeatMasker installation folder:**
- DupMasker
- ProcessRepeats
- RepeatMasker
- RepeatProteinMask

**Scripts in RepeatMasker/util folder:**
- buildRMLibFromEMBL.pl
- buildSummary.pl
- calcDivergenceFromAlign.pl
- createRepeatLandscape.pl
- getRepeatMaskerBatch.pl
- queryRepeatDatabase.pl
- queryTaxonomyDatabase.pl
- RM2Bed.py
- rmOut2Fasta.pl
- rmOutToGFF3.pl
- rmToUCSCTables.pl
- trfMask
- wublastToCrossmatch.pl
