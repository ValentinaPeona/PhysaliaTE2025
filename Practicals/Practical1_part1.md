# Practical 1 - part 1: De-novo repeat prediction with RepeatModeler2

*Valentina Peona*

*valentina.peona@nrm.se*

*valentina.peona@gmail.com*

## Introduction

Obtaining a high quality repeat library that encompasses the entire diversity of repeats present in the genome of interest is key to get reliable downstream analyses related to transposable elements.

If you work on model organisms, it is  likely that the diversity of transposable elements of your species has already been characterised and that you can find a library of consensus sequences avaialable on databases like (Repbase)[https://www.girinst.org/repbase/] and (Dfam)[https://www.dfam.org]. In this case, it is not essential to run tools like `RepeatModeler2`. This can be true also if your species of interest is closely related to a species whose transposable elements have been already characterised. However, we find over and over that masking a genome of a new species using repeat libraries from other species can be insufficient. This happens when the genome of interest contains species-specific repetitive elements never described before and these elements would remain largely unmasked or partially and/or incorrectly annotated. Furthermore, curated repeat libraries from sister species can actually reciprocally help the annotation of the two genomes [Boman et al. 2018](https://www.mdpi.com/2073-4425/10/4/301).



📝 Libraries and consensus sequences 📝

To get a complete overview of the diversity of transposable elements in the genome of interest then we need to run a **de-novo repeat prediction** that will produce a new **repeat library made of raw consensus sequences**. These raw consensus sequences are going to need a good round of manual curation (`Practical 2`) and after that the library will be ready to annotate our genome of interest.

> A TE family (that can be seen as a species) can be represented by a consensus sequence approximating that of the ancestral progenitor. Such consensus sequence can be created on top of a multiple alignment of individual genomic copies (or “seeds”) from which each ancestral nucleotide can be inferred based on a majority rule along the alignment. Similarly, the seed alignment may be used to generate a profile Hidden Markov Model (HMM) for each family. Flynn et al., 2020
A library of consensus sequences will then be essential to find sequences in the genome of interest that are similar to the consensus sequences themselves, namely the different copies of the TEs.

A good repeat library has these characteristics:

🟢 is complete - the entire diversity of repeats is represented

🟢 contains nonredundant consensus sequences - each element is represented only once

🟢 contains full-length consensus sequences - each elements is not fragmented/truncated

There are a few tools available for the de-novo characterisation of repeats (e.g., `CARP`, `REPET`, `EarlGrey`, `RepeatExplorer2` and other all listed in the [TE Hub website](https://tehub.org)) but for time reasons we will use only `RepeatModeler2`.

To run `RepeatModeler2` we need a genome of interest (target genome) and we will use the assembled genome of the tardigrade *Ramazzottius varieornatus* (GCA_001949185.1; abbreviated throughtout the practicals as `ramVar`; [Yoshida et al. 2017](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2002266)). 

<p align="center"> 📖 <strong>Learning objectives</strong> 📖 </p>

- [ ] create of a de-novo repeat library for `ramVar` using `RepeatModeler2`
- [ ] inspect the outout of `RepeatModeler2` and describe what a repeat library and consensus sequences are

#### ⚠️ Note
Given the short time within the course, you will not be able to fully run the `RepeatModeler2` commands, therefore we generated all the output for you that can be found in the folder `~/Share/TE25/Practical1/`. We suggest you to follow the tutorial and try the commands out but kill the analysis right after. All the output necessary is available to you in any moment!


---
# Part 1

## Getting ready

:books: **Important bash commands**

You don't need to be an expert of Unix and Bash to follow this tutorial, however here I give you a short list of essential bash commands that can be useful for the practicals. 

```bash
# change directory
cd
# get the name of the folder in which you are located
pwd
# list files in directory
ls
ls -l
# visualise a file
less <filename>
# exit less by typing q
# print the content of a file
cat <filename>
# print first lines of a file
head <filename>
# print last lines of a file
tail <filename>
# edit a file
nano <filename>
# create directory
mkdir <dir_name>
# copy file
cp /path/to/file/file.txt ./
# create symbolic link
ln -s /path/to/file/file.txt ./
# get help for how to run commands
man <command>
# if man does not work, try with
<command> -h
```

:open_file_folder: **Create your working directory**

Create your own folders

```bash
# create folders
mkdir -p ~/TE25/Practical1/Data ~/TE25/Practical1/Code ~/TE25/Practical1/RMDL

# copy reference genome and repeat library in the Data folder
cp ~/Share/TE25/Practical1/Data/* ~/TE25/Practical1/Data/
cp ~/Share/TE25/Practical1/Code/* ~/TE25/Practical1/Code/
```

---

## Run RepeatModeler2

A run of `RepeatModeler2` on a genome consists of two steps:
1. indexing of the genome (somewhat similar to what bwa index does) using the command `BuildDatabase`
2. clustering of similar sequences to create the de-novo library using the command `RepeatModeler`

<p align="center"><figure><img src="./repeatmodeler2.jpg" alt="RepeatModeler2 pipeline" /><figcaption><strong>Figure 1.</strong> RepeatModeler2 pipeline. Image from <a href="https://www.pnas.org/doi/10.1073/pnas.1921046117">Flynn et. al 2020</a>.</figcaption></figure></p>

The core of the `RepeatModeler2` pipeline (illustrated above in **Figure 1**) consists of up to 6 rounds of sequence subsampling of your genome. Each round subsamples chunks of your genome of different sizes and on each of these chunks, several tools are run to identify both tandem and interspersed repeats. At the end of each round, multisequence alignments of repetitive elements are created and consensus sequences are produced on top of them and stored in the `consensi.fa` file. `RepeatModeler2` can also run an optional additional pipeline called `LTRstruct` particularly dedicated to improve the characterisation of LTR retrotransposons.

On the Amazon server we installed `RepeatModeler2` in a conda environment that you can activate like this:

```bash
conda activate te25
singularity exec ~/Share/TE25/dfam-tetools-latest.sif RepeatModeler
```

Let's give it a try!

```bash
# go to working directory
cd ~/TE25/Practical1
```

As first thing call the command `BuildDatabase`, look at the options and run it on `ramVar.fasta` genome.

<details>
    <summary><strong>BuildDatabase command</strong></summary>

```bash
# STEP1: create a database for the genome. We put the database in the Data folder together with the genome assembly
# BuildDatabase takes about 1 minute to run
# BuildDatabase [-options] -name <name of the database> <genome file>
BuildDatabase -name Data/ramVar Data/ramVar.fasta
```

</details>

Now that your index/database is ready, we can run `RepeatModeler2` on it.

Call the `RepeatModeler` command and inspect the options, do you see anything useful/interesting?
- What input files do you need to run RepeatModeler2?
- Which options would you use?
- What does `LTRStruct` mean/do?

Also go to the method section of the paper of [RepeatModeler2](https://www.pnas.org/doi/full/10.1073/pnas.1921046117#sec-1) and read the about the LTR module of the tool.

Run RepeatModeler2 on the tardigrade genome.

<details>
    <summary><strong>RepeatModeler2 command</strong></summary>

```bash
# STEP2: run RepeatModeler2 on the database you just created
RepeatModeler -database Data/ramVar
```

⚠️ The conda installation on the AWS doesn't include the software for the LTRStruct!

</details>

:warning: **Warning**

Once you managed to start `RepeatModeler2`, let it run for a minute or two (or until you lose your patience) to see the very first round of repeat prediction start and then please kill the command with `Ctrl + C` as we do not have the time to let it finish. You will see that a folder with a name similar to this `RM_1962636.MonJan131517222025` has been created in your working directory that contains the intermediate files of the RepeatModeler2 pipeline. You can give it a look and then you can jump directly to the final output files of the pipeline.

:exclamation: Copy the output from the `~/Share/TE25/Practical1/RMDL` folder to yours:

As you can see there are two commands to run to get the final output from the tool (the library of raw consensus sequences). The first command creates a database of your genome assembly fasta file, basically it indexes the fasta to better access it during the de-novo characterisation similarly to what BLAST does. Indeed, most of the `RepeatModeler2` analysis consists of a large number of alignments.

#### 📌 Note
When you're going to run `RepeatModeler2` and `RepeatMasker` on your fasta files, pay attention to the presence of strange/special characters in your fasta headers, they may cause problems during the analysis. Also, it is suggested to have rather short fasta headers or they will be cut and you may lose some important pieces of information: this is of particular importance when running RepeatMasker.

❗**Copy the output files in your own RMDL folder**
```bash
cp -r ~/Share/TE25/Practical1/RMDL ~/TE25/Practical1/RMDL
```

---

## Inspect the RepeatModeler2 output

📖 **Let's give a look at the output of RepeatModeler2** 📖

Go to your `RMDL` folder and list the files present there.

:question::question: **Questions** :question::question:

🔴 What files do you find? What is the difference between a fasta file and a Stockholm file?

🔴 What is the difference between the classified and non-classified fasta files?

🔴 What is a consensus sequence?

🔴 Do you see shortcomings in the use of consensus sequences to annotate repeats?

🔴 Do you think this library is complete, non-redundant and contains full-length elements? (Just reflect on this, it is difficult to know at this stage.)

🔴 What types of repeats are found in your library? Use the classified version of the library to answer this.

🔴 Do you think this library contains only transposable elements?

Run this `librarySummary. R script to get more info about the characteristics of the library. To run it, the `consensi.fa.classified` file must be converted from multiline fasta to single line fasta.

```bash
cd ~/TE25/Practical1/RMDL

# convert fasta to single line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < consensi.fa.classified | awk 'NR > 1' > ramVar_rm1.0_temp.fasta

# rename the headers: replace "rnd" with a custom name, in this case "ramVar"
perl ~/TE25/Practical1/Code/renameRMDLconsensi.pl ramVar_rm1.0_temp.fasta ramVar ramVar_rm1.0.fasta
rm ramVar_rm1.0_temp.fasta

Rscript ~/TE25/Practical1/Code/librarySummary.R ramVar_rm1.0.fasta
```

<details>
    <summary><strong>Answers</strong></summary>

<strong>What files do you find?</strong>
In the <code>Results</code> folder you find the final output of RepeatModeler2. Note that the intermediate files have been removed. RepeatModeler2 creates a huge amount of intermediate files that are useful essentially to resume a run that failed.
<br>
The fasta and Stockholm files contain the consensus sequences generated by the tool. While the fasta file can be directly used by RepeatMasker, the Stockholm file can be used to create HMM profiles and be used by other HMM-based tools.
<br>

```
consensi.fa  : Consensus sequences
consensi.fa.classified: Consensus sequences classified
families.stk : Seed alignments
families-classified.stk: Seed alignments classified
```

<br>

<strong>What is the difference between the classified and non-classified fasta files?</strong>
The difference lies in the headers of the fasta files. In the classified version you find a tag after the name of the repeat that specifies (when possible) the type of repeat. Note that this classification can be partial and sometimes incorrect.
<br><br>

<strong>What is a consensus sequence? Do you see shortcomings in the use of consensus sequences to annotate repeats?</strong>
A consensus sequence is an approximation of the ancestral repeat that gave rise to the insertions found in the genome. Note that a consensus sequence is made on top of a set of similar (and likely homologous) repetitive elements, therefore it represents an approximation of the ancestor of that limited set of aligned sequences therefore could not capture the entire variability of that family of TE.
<br><br>

<strong>Do you think this library is complete, non-redundant and contains full-length elements?</strong>
Most likely, the repeat library just produced here does not meet all these criteria. Likely it is incomplete because low copy number repeats are still difficult to identify and build a consensus on top of a very few sequences. Likely some of the consensus sequences are not full-length because of the nature of the transposable elements (e.g., LINEs and 5' truncation) or because of the fragmentation of the genome used. However, the library should rather be non-redundant because RepeatModeler2 has a step in the pipeline that gets rid of most of the redundancy. It is good practice to check for it! <br> A library can be improved with some manual curation!
<br><br>

<strong>Do you think this library contains only transposable elements?</strong> The library likely contains also tandem repeats but it can also contain multi-copy genes. RepeatModeler2 takes care of removing gene-related sequences but it is still possible that some genes remain. You can check for it using blast or tools like <a href="https://github.com/NBISweden/ProtExcluder">ProtExcluder</a>.


</details>

:note: Note that the header names of the consensus sequences in the `consensi.fa.classified` follows a precise nomenclature and style which is necessary for RepeatMasker to give you correct abundance estimates for each kind of TE.

---

<p align="center"> 🔎 <strong>Further readings and tools</strong> 🔎 </p>

- [A beginner’s guide to manual curation of transposable elements](https://link.springer.com/article/10.1186/s13100-021-00259-7). This guide includes several protocols, videos and code.
- [Earl Grey: A Fully Automated User-Friendly Transposable Element Annotation and Analysis Pipeline](https://academic.oup.com/mbe/article/41/4/msae068/7635926). Example of a tool to get an automatically curated repeat library.

---
