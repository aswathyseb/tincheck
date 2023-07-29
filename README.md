# Tincheck

Tincheck is a python package to calculate transcript Integrity Number (TIN) and Transcription overlap.

Transcripts with uneven coverage or transcription overlap can result in false positives in differential expression analysis.
Transcript Integrity Number is a metric that calculates coverage evenness and can be used as a filtering criteria in RNA-Seq studies to improve the accuracy of results.

### Quick Start

Test data and annotation file is in `data` folder.

**How to run**
    
    tincheck tin -a data/ann.gtf data/sample.bam
    
**Inputs**

1. Alignment file in bam format
2. Annotation file in GTF/GFF3 format

**Output**

Tab delimited textfile  with TIN score calculated for each gene/transcript/any other feature specified in the input.
An example is given below.

    target_id       eff_length  S1_count S1_exp_tin  S1_obs_tin
    PF3D7_0102700	1683	    670	        100.0	    70.9
    PF3D7_0103700	1624	    135	        100.0	    72.8
    PF3D7_0107300	1581	    4508        100.0	    70.4
    PF3D7_0107600	5702	    4979        100.0	    74.9
    PF3D7_0107800	4424	    924	        100.0	    78.8


**How to install the script?**

Step1: Install script

    pip install tincheck --upgrade

Step2 :Install additional requirements
    
    conda install --file conda-requirements.txt

### Additional details


TIN denotes how evenly a feature is covered by reads. By default, the script calculate the coverage evenness across a 
gene by considering the coverage across all exons. However, the script can calculate coverage evenness across a transcript or any other feature that is in the annotation file.
The annotation file should have a gene feature row present in it.

**How to calculate tin for each gene?**

    tincheck tin --a data/ann.gtf data/sample.bam 

**How to calculate tin for each transcript?**

    tincheck tin -g transcript_id -a data/ann.gtf data/sample.bam 

Here -g option should be the transcript grouping attribute present in the annotation file. 


**How to calculate TIN across coding regions of a gene?**
    
    tincheck tin -f CDS -a data/ann.gtf data/sample.bam


**How to calculate TIN across coding regions of a transcript?**
    
    tincheck tin -f CDS -g transcript_id -a data/ann.gtf data/sample.bam


**How is tin calculated**

Overlapping features specified by `--feat` are merged together and coverage evenness is captured using Shannon's entropy formula  (H).

This is then converted into TIN score as

    TIN = (100*exp(H))⁄length  where H= Shannon's entropy formula.

**Transcription overlaps**

Transcription overlaps are flagged using `overlap.py` script.

	tincheck overlap --ann data/ann.gtf data/sample.bam >sample_overlap.txt

By default, genes are checked for overlap only if 

    gene-tin < tin-cutoff and gene-count > count-cutoff
    
 ie, genes with enough counts but low tins are checked for transcription overlap
 
A gene is flagged to have a transcription overlap if either of its neighboring genes and the corresponding inter-genic region have a higher tin-score and read count than the gene of interest.

In stranded mode, if the count and tin values of neighboring genes and the intergenic region in the gene sense strand is 
comparable to the gene count and tin, then that gene is considered to have an overlap from neighboring transcripts.


More specifically in stranded mode, a gene is considered to have a neighboring transcript overlap, if **either of** the following conditions are satisfied.

    1. left gene tin in gene sense strand >= gene tin and left intergenic tin in gene sense strand >=gene-tin
    2. right gene tin in gene sense strand >= gene tin and right intergenic tin in gene sense strand >=gene-tin

