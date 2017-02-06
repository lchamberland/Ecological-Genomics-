

## Ecological Genomics

## *2016.02.06* - RNAseq

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



# *2016.02.06* - coding - RNAseq data  

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

to run file 

```
./filename or bash.filename
```

Output for my data files 

*[lchambe1@pbio381 scripts]$ ll*

*total 4*

*-rwxr--r--. 1 lchambe1 users 809 Feb  6 11:33 trim_example.sh*

[lchambe1@pbio381 scripts]$ ./trim_example.sh 

*TrimmomaticPE: Started with arguments: -threads 1 -phred33 /data/project_data/fastq/105-11H0R1.fq.gz /data/project_data/fastq/105-11H0R2.fq.gz /data/project_data/fastq/cleanreads/105-11H0R1clean_paired.fq /data/project_data/fastq/cleanreads/105-11H0R1clean_unpaired.fq /data/project_data/fastq/cleanreads/105-11H0R2clean_paired.fq /data/project_data/fastq/cleanreads/105-11H0R2clean_unpaired.fq ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:6:28 HEADCROP:9 MINLEN:35*

*Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'*

*ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences*

*Input Read Pairs: 23373711 Both Surviving: 19911093 (85.19%) Forward Only Surviving: 2539010 (10.86%) Reverse Only Surviving: 356754 (1.53%) Dropped: 566854 (2.43%)*

*TrimmomaticPE: Completed successfully*

*[lchambe1@pbio381 scripts]$* 