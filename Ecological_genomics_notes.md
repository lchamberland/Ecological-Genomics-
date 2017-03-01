# Online Notebook

### Author: Lisa Chamberland

### Ecological Genomics:

## Overall Description of notebook

fill in your description here!

## Date started: (2017-02-06)

## Date end:   (2017-02-15)

## Philosophy

Science should be reproducible and one of the best ways to achieve this is by logging research activities in a notebook. Because science/biology has increasingly become computational, it is easier to document computational projects in an electronic form, which can be shared online through Github.    

### Table of contents for 10 entries (Format is *Page: Date(with year-month-day). Title*)

- [Page 1: 2017-02-06](#id-section1). RNAseq

- [Page 2: 2017-02-08](#id-section2). Transcriptomics

- [Page 3: 2017-02-13](#id-section3). Transcriptomics 2

- [Page 4: 2017-02-15](#id-section4) Transcriptomics 3

- [Page 5:](#id-section5) Transcriptomics 4

- [Page 6:](#id-section6) 

- [Page 7:](#id-section7).

- [Page 8:](#id-section8).

- [Page 9:](#id-section9).

  ___________

  <div id='id-section1'/>

## Page 1: 2017-02-06. RNAseq

### Info update - Melissa 

#### *De Wit et al., 2012* 

*2016.02.06*

### RNAseq

* Approach 
* Experimental design 
* Library prep 
* Sequencing (facility)
* Recieve the data 
* Computer/server setup 

**1. Clean reads** (evaluate quality) *fastq files*

* adapters
* nucleotide quality 
* length 

**2. Evalutate Quality**

**3. de novo transcriptome assembly** *fasta*

* (ref)
* evaluate your assembly 
  * compare to closely related species 
  * compare to a reference set of genes (Core Eukariotic Genes (KEG/CEG))
  * N50
  * number of contigs 
* Annotation (BLAST search)

**4. Map reads to reference transcriptome** *SAM - sequence alignment file* 

* alignment files 

**5a. Extract read count info**

* number of reads that map to each contig for each sample

**5b. Identify SNPs**

**6a. DGE analyzes**

* co-exp network analyses

**6b. Population genetics**

* genetic difference 
* population structure 



### Paper Discussion - Muhammhad 

#### *Dunning et al., 2014*

50bp single end short reads for their RNA seq

​	harder to assemble if you are working from scratch, but you get more of them 

* assemble de novo because no reeference transciptome 
* major analyses
  * differential expression of RNA seq data 
* **Trinity** - collapse splice varient contigs into "unigenes" - assembly stage 
  * it will take the longest gene out of the cluster of the alternatively spliced genes
* Differential expression 
  1. DESeq
  2. EdgeR
  3. BaySeq
* Enrichment analyis
  * test the GO (Gene Ontology) terms associated with differentially regulated cold-responsive unigenes and those associated with non-differentially expressed unigenes
    * Group genes into pathways to make sense of the response to the particular treatments 
    * *Is this something Hannah and I can do for our project if we get many genes back?*
    * BLASTX - Translate genes to proteins and blast 
      * nr - gene annotation 
      * uniprot database - Gene Ontology (GO) functional enrichment 
        * Given what is in the functional categories in the transcriptone are my genes non-random? 
        * Candidate genes - came out in the test as interesting 
          1. Differentially expressed 
          2. part of group that are over expressed 
  * Biological replicates - three different individuals 
    * how consistent these are among your treatment 
  * Technical replicates - form individual 1 you take 3 samples from one individual - expect these to very precise 
  * most of the up and down regulation seem to be unique to a given focal population 
  * qPCR (quantitative PCR) - looking at copy number 
* Are these expression differences driven by genetic drift or selection?
* Choose statistal test and stick with it or make it a comparison 


### Coding 

login to server 

```
ssh lchambe1@pbio381.uvm.edu
```

enter password

enter data folder

``` 
cd /data
```

list files 

```
ll
```

enter subfolder 

```
cd project_data/
```

My file(s) that I will be cleaning and saving to a common directory: 

* 10 5-11 H 0 R1.fq.gz
* 10 5-11 H 0 R2.fq.gz

Using the program **zcat** to look at the zipped file 

```
zcat filename | head 
```

Lines 

1. A unique identifier about sequencer used, it's cluster, adapter etc 
2. Sequence data
3. +
4. Sequence of quality scores for each base - represents probability of error (Q score) 
   * probability that nucleotide is erronious 
   * Higher q score better 
   * quality score greater than or equal to 30 



Use program **fastqc** to clean the data file 

copy example file into your home file and folder scripts 

```
cp trim_example.sh ~/scripts/
```

peak at first couple lines

```
head filename 
```

manuplate file using javascript file

```
vim filename
```

```
java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
                -threads 1 \
                -phred33 \
                 /data/project_data/fastq/10_5-11_H_0_R1.fq.gz \
                 /data/project_data/fastq/10_5-11_H_0_R2.fq.gz \
                 /data/project_data/fastq/cleanreads/"10_5-11_H_0_R1_clean_paired.fq" \
                 /data/project_data/fastq/cleanreads/"10_5-11_H_0_R1_clean_unpaired.fq" \
                 /data/project_data/fastq/cleanreads/"10_5-11_H_0_R2_clean_paired.fq" \
                 /data/project_data/fastq/cleanreads/"10_5-11_H_0_R2_clean_unpaired.fq" \
                 ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
                 LEADING:28 \
             TRAILING:28 \
             SLIDINGWINDOW:6:28 \
             HEADCROP:9 \
             MINLEN:35 \
~                                             
```

Above is an example of how you need to edit this example file. 

* You are adding your R1 and R2 filenames

* change .fa to .fq

* change path name 

* filenames and path

* input (2)

* output (4) - also .fq

  ​

To save a file and quit

```
escape
:w
:q
```

file is executable if it is lit up green

to run/execute file 

```
./filename or bash.filename
```

Output for my data files 

*[lchambe1@pbio381 scripts]$ ll*

Updated upstream

*total 4*

*-rwxr--r--. 1 lchambe1 users 809 Feb  6 11:33 trim_example.sh*

[lchambe1@pbio381 scripts]$ ./trim_example.sh 

*TrimmomaticPE: Started with arguments: -threads 1 -phred33 /data/project_data/fastq/105-11H0R1.fq.gz /data/project_data/fastq/105-11H0R2.fq.gz /data/project_data/fastq/cleanreads/105-11H0R1clean_paired.fq /data/project_data/fastq/cleanreads/105-11H0R1clean_unpaired.fq /data/project_data/fastq/cleanreads/105-11H0R2clean_paired.fq /data/project_data/fastq/cleanreads/105-11H0R2clean_unpaired.fq ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:6:28 HEADCROP:9 MINLEN:35*

*Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'*

*ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences*

*Input Read Pairs: 23373711 Both Surviving: 19911093 (85.19%) Forward Only Surviving: 2539010 (10.86%) Reverse Only Surviving: 356754 (1.53%) Dropped: 566854 (2.43%)*

*TrimmomaticPE: Completed successfully*

*[lchambe1@pbio381 scripts]$* 



scp 

______

<div id='id-section2'/>

## Page 2: 2017-02-08. Transcriptomics

### Paper Discussion - Laura 

De Panis, D.N., Padró, J., Furió‐Tarí, P., Tarazona, S., Milla Carmona, P.S., Soto, I.M., Dopazo, H., Conesa, A. and Hasson, E., 2016. Transcriptome modulation during host shift is driven by secondary metabolites in desert Drosophila. *Molecular Ecology*, *25*(18), pp.4534-4550.



Experimental design

* samples collected in environment
* biological replicates - 3 lines - 3 different inversions 
* 4 treatments - 2 assess host plant effect, 2 assess



### Coding 

Goals today

* finish cleaning (trimmomatic)
* fastqc(visualize)
* make table of number of reads
* design assembly tests 
* start assemblies
* evaluate assembly 
* Script 
  * paths 
    * program 
    * input 
    * output
  * filename
    * in 
    * out 
* Moving through directories
  * moving files 
  * scp - secure copy to move from server to your pc 
* executing scripts 
* calling program 



to rename old file name to a new filename 

```
mv old_filename new_filename 
```

to clean and create output html file 

* *this will put your cleaned file in whatever directory you are in*

```
fastqc filename.fq.gz
```

Move .html file to your computer desktop using the **scp** command (*open another terminal*)

```
scp lchambe1@pbio381.uvm.edu:/data/project_data/fastq/filename.html ~/Desktop/
```

What will affect making our assembly 

* using paired reads 
* number of reads 

Starting a de novo assembly using Trinity 

```
$ /data/program/trini
```
>>>>>>> Stashed changes

*total 4*

*-rwxr--r--. 1 lchambe1 users 809 Feb  6 11:33 trim_example.sh*

[lchambe1@pbio381 scripts]$ ./trim_example.sh 

*TrimmomaticPE: Started with arguments: -threads 1 -phred33 /data/project_data/fastq/105-11H0R1.fq.gz /data/project_data/fastq/105-11H0R2.fq.gz /data/project_data/fastq/cleanreads/105-11H0R1clean_paired.fq /data/project_data/fastq/cleanreads/105-11H0R1clean_unpaired.fq /data/project_data/fastq/cleanreads/105-11H0R2clean_paired.fq /data/project_data/fastq/cleanreads/105-11H0R2clean_unpaired.fq ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:6:28 HEADCROP:9 MINLEN:35*

*Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'*

*ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences*

*Input Read Pairs: 23373711 Both Surviving: 19911093 (85.19%) Forward Only Surviving: 2539010 (10.86%) Reverse Only Surviving: 356754 (1.53%) Dropped: 566854 (2.43%)*

*TrimmomaticPE: Completed successfully*

*[lchambe1@pbio381 scripts]$* 

______

<div id='id-section3'/>

## Page 3: 2017-02-13. Transcriptomics 3 

### Info Update - Lauren

**Glossary**

- sequence coverage - the average nuber of reads that align/"cover" known reference bases
- read depth - total numer of bases sequenced/aligned at given reference base position 
- transitional noise - unexplained variation/ randomness 
- power - probability of rejecting false null hypothesis
- biological variation - natural variation in gene 

Background

* enables DE examination (inter-population individual)

  * disease resistance
  * mating behavior 
  * adaptive significance

* molecular mechanisms —> phenotypic/ behavioral plasticity, migration patterns

* limitations 

  * ref. genome quality 
  * gene annotation availability
  * expence per sample 
  * lib preb

* Issues 

  * under utilization of biological replicates
    * requiring indep library preparations
    * doesn't include pooled samples
    * 23/158 studies (%15) > 3 biological reps
    * devise broad bio conclusions
  * prioritize seq depth over replication <— problem 
  * wide dynamic range of RNA - seq data 
    * noisy 
      * poisson counting error
      * technicla variance
      * biological variance 
        * **lower variance gives you higher power** 

  General rules of thumb

  1. use more bio. replicates instead of depth

  2. sequence depth > 10 reads/ transcript 

     a. ~10 - 20 M mapped reads/samples

  3. 3 biological replicates per conditions 

  4. conduct a pilot experiment 

     a. what is best/ powerful experiment I can afford?

     b. what is the smallest fold change I can detect

  ### Paper Discussion - Sam 

  *Johnston et al. 2016 Seasonal gene expression in migratory songbird* 

* Tissue matters 

* statistical issues

  * number of reps
  * power coming from the 31-48 million reads 

### Coding

Open reading frame - starts with a start codon and ends with a stop codon 

* used to find complete transcripts and the longest ones

bwa - program for mapping reads 

quit/stop a program

```
ctrl c 
```





```
screen 
```

 ```
screen - r 
 ```

______

<div id='id-section4'/>

## Page 4: 2017-02-15. Transcriptomics 4

### Info Update - Sam

SNPs and population genomics 

SNP data - expressed sequences

tissue -> sequence -> clean/trim -> assembly -> SNP detection/ validation => practical applications

1. Tissue
   * breadth of tissue, developmental stages because of exon skipping
2. pool and sequence libraries
   * ~30 -100 m paired end long reads
3. process raw sequence data
   * important for the SNP detection 
4. digital normilization 
   * remove high coverage reads and associated errors 
   * loss of quantitative info 
5. Assemble cleaned pair end long reads
6. Prune 
   * reduce DNA contamination, non coding RNA, gene fragments 
7. Assembly evaluation 
   * reference genome if you have on e
   * COGS - concerved eukaryotic ortholoc genes 

SNP detection

* software - constant patters of sequenc variation 
  * Problems
    * sequence error - software eliminates SNPs of low frequency 
    * artifacts caused by INDELs (insertions or deletions)
      * filter SNP clusters near INDELS 
      * quality scores 
* Validation 
  * primers 
  * sequencing using mas spec 

Applications

* differences in population structure
* natural selection acting on a particular loci 

1. outliers - for a given locus, what's the level of differentiation compared to differentiation accross the genome - looking at Fst values that are at the end of the normal distribution and assumes directional selection 
2. non-outlier - tests high Fst loci for other features associated with selection 
   * fitness 
   * enrichment for in functional roles 

RPKM - reads per kilobase transcript per million

* normalizing - scaling and comparing them based on expression instead of size 



References for samples that I need to check.

```
10_5-11_H_0_R1.cl.pd.fq  10_5-11_H_0_R1.cl.un.fq  
```

### Coding

1. Lab notebook 
2. SAM files
3. Extract expression 
4. What's going on in the background 



-rw-r--r--. 1 lchambe1 users  4887950297 Feb 15 10:21 10_5-11_H_0_R1_clean_paired.fq

-rw-r--r--. 1 lchambe1 users   571278338 Feb 15 10:21 10_5-11_H_0_R1_clean_unpaired.fq

-rw-r--r--. 1 mpespeni users  5247438419 Feb 10 11:15 10_5-11_H_0_R1.cl.pd.fq

-rw-r--r--. 1 mpespeni users   363350299 Feb 10 11:15 10_5-11_H_0_R1.cl.un.fq

-rw-r--r--. 1 lchambe1 users  4881390207 Feb 15 10:21 10_5-11_H_0_R2_clean_paired.fq

-rw-r--r--. 1 lchambe1 users    76576027 Feb 15 10:21 10_5-11_H_0_R2_clean_unpaired.fq



s=search, find ::, replace with _, g is global - replace globally 

```
sed -i 's/::/|_/g' YOURFILE.sam
```



error

```
[bwa_seq_open] fail to open file '/data/project_data/fastq/cleanreads/10_5-11_H_0_R1.fq.gz_left_clean_paired.fq' : No such file or directory
```



```
#!/bin/bash 

# To run from present directory and save output: ./bwaaln.sh > output.bwaaln.txt 

myLeft='10_5-11_H_0_R1.fq.gz_left_clean_paired.fq'
echo $myLeft

myRight=${myLeft/_R1.fq.gz_left/_R2.fq.gz_right}
echo $myRight

myShort=`echo $myLeft | cut -c1-11`
echo $myShort

# bwa index /data/project_data/assembly/longest_orfs.cds  # This only needs to be done once on the reference

bwa aln /data/project_data/assembly/longest_orfs.cds /data/project_data/fastq/cleanreads/$myLeft > $myLeft".sai"
bwa aln /data/project_data/assembly/longest_orfs.cds /data/project_data/fastq/cleanreads/$myRight > $myRight".sai"
bwa sampe -r '@RG\tID:'"$myShort"'\tSM:'"$myShort"'\tPL:Illumina' \
        -P /data/project_data/assembly/longest_orfs.cds $myLeft".sai" $myRight".sai" \
        /data/project_data/fastq/cleanreads/$myLeft \
        /data/project_data/fastq/cleanreads/$myRight > $myShort"_bwaaln.sam"
~                                                                              
```



<div id='id-section5'/>

## Page 5: 2017-02-17. 

/data/project_data/assembly/longest_orfs.cds



<div id='id-section6'/>

## Page 6: 2017-02-27. Guest Speaker (Scott Edwards)

***Glossary***

**coalescent** - common ancestor- node - when you sequence many individuals and you find a pattern of where they are coalescing, you can make an inference 

**reticulation** describes the origination of a [lineage](https://en.wikipedia.org/wiki/Lineage_(evolution)) through the partial merging of two ancestor lineages, leading to relationships better described by a [phylogenetic network](https://en.wikipedia.org/wiki/Phylogenetic_network) than a bifurcating [tree](https://en.wikipedia.org/wiki/Phylogenetic_tree).

**purifying/background selection**

**gene trees vs. species**

**introgression/recombination** - requires a history of divergenc; gene flow; hybridization 

**incomplete lineage sorting (ILS)** – so®†ing of taxa is no† comple†e ye†, 

* large population sizes and short time will result in genes that are still around
* time is hight of the tree, population size is the width of the brach 
* when Ne is very big, ratio is very small.  Ne is population size 
* t/Ne
* Genomes are collections of evolutionary histories





```
allcountsdata.txt                             100% 3883KB   3.8MB/s   00:00    
cols_data_trim.txt                            100% 2342     2.3KB/s   00:00    
countsdata_trim2.txt                          100% 3394KB   3.3MB/s   00:00    
countstatsummary.txt                          100% 8059     7.9KB/s   00:00    
DESeq2_SSW_round2.R                           100%   13KB  13.3KB/s   00:00    
scp: /data/project_data/DGE/round1: not a regular file
[lchambe1@pbio381 Desktop]$ cd ~/Desktop/sea_star2
```



<div id='id-section5'/>

## Page 7: 2017-03-01 

**DESeq2** 

LRT - Likelihood ratio test - comparing two different models to see hwich one is better - does this model explain the significance in the data and how accurate is it 

* log(likelihood(full model) - likelihood(reduced model)) = Likelyhood ratio (LR)
* Health, day, health*day - health, day 
* distributed as chi square with df = # parameters full model - parameters reduced
* given the data this is the most likely model 
* less parameters = more power ; simplest model 

Model 3 

*give group variable that you have defined as all combinations of your variables.  For each of the individuals*

PCA plot - clusters on arbitrary axis based on similarites etween two samples 

*you cant have individuals and day in the same model because you days are the replicates in the model*

Moving files from the server to your desktop 

```
scp lchambe1/data/project_data/DGE/*.csv . 
```

```
scp lchambe1/data/project_data/DGE/*.txt . 
```

