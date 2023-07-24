# Tincheck

Tincheck is a python package to calculate transcript Integrity Number (TIN) and Transcription overlap.

Transcripts with transcription overlap or uneven coverage can result in false positives in differential expression analysis.
Transcript Integrity Number is a metric that calculates coverage evenness and can be used as a filtering criteria in RNA-Seq studies to improve the accuracy of results.


**Transcript Integrity Number (TIN)** calculation

TIN denotes how evenly a feature is covered by reads. By default, the script calculate the coverage evenness across a 
gene by considering the coverage across all exons. However, the script can calculate coverage evenness across a transcript or any other feature that is in the annotation file.


**Inputs**
1. Alignment file in bam format
2. Annotation file in GTFor GFF3 format

**Output**

Tab delimited textfile  with TIN score calculated for each gene/transcript/any other feature specified in the input.
An example is given below.

    target_id       eff_length  S1_count    S1_exp_tin  S1_obs_tin
    PF3D7_0102700	1683	    670	        100.0	    70.9
    PF3D7_0103700	1624	    135	        100.0	    72.8
    PF3D7_0107300	1581	    4508        100.0	    70.4
    PF3D7_0107600	5702	    4979        100.0	    74.9
    PF3D7_0107800	4424	    924	        100.0	    78.8


**How to install the script?**
    
    python setup.py install

It would be better to create a conda environment and install the script as shown below
    
    conda create -n tin -y python=3.6
    conda activate tin

Install additional requirements
    
    conda install --file conda-requirements.txt 

Finally, install the script in the conda environment as
    
    python setup.py install

**How to run the script?**

The required inputs to the script is bam file(s) and annotation file in GTF or GFF format.
 
The annotation file should have a gene feature row present in it.

Given a bam file and an annotation file, the simplest way to run the script is

    tincheck tin --ann data/ann.gtf data/sample.bam 

**How to calculate TIN across coding regions of a gene?**
    
    tincheck tin  --ann data/ann.gtf --feat CDS data/sample.bam 
 

**Available options**

`--ann` Annotation file GTF/GFF3 format.

`--feat` Subfeature that is grouped together to compute TIN. It should match the 3rd column in annotation file.
         Default is exon.
         
`--groupby` Feature by which subfeatures needs to be combined. It should match the 3rd column in annotation file.
            Default is gene.
         
`--strand` When specified, depth is calculated in strand-specific manner. Possible values are 'both', 'sense' or 'antisense'. Default is both.


`--lib_type` Specify if the library is paired or single. Valid options are 'single' or 'paired. Default is single.'


`--bg` If specified background noise will be subtracted
 
`--n` No. of bases to be subtracted from each ends of the feature to calculate effective length. Default is 50.
    
**How is tin calculated**

Overlapping features specified by `--feat` are merged together and coverage evenness is captured using Shannon's entropy formula  (H).
This is then convered into TIN Score as

    TIN = (100*exp(H))⁄length  where H= Shannon's entropy formula.


**Transcription overlaps**

Transcription overlaps are flagged using `overlap.py` script.

	tincheck overlap --ann data/ann.gtf data/sample.bam >sample_overlap.txt

By default, genes are checked for runins only if 

    gene-tin < tin-cutoff and gene-count > count-cutoff
    
 ie, genes with enough counts but low tins are checked for transcription overlap
 
If the count and tin values of neighboring genes and the inter-genic regions in the gene sense strand is 
comparable to the gene count and tin, then that gene is considered to have an overlap from neighboring transcripts.


More specifically, a gene is considered to have a neighboring transcript overlap, if **either of** the following conditions are satisfied.

    
    1. left gene tin in gene sense strand >= gene tin and left intergenic tin in gene sense strand >=gene-tin
    2. right gene tin in gene sense strand >= gene tin and right intergenic tin in gene sense strand >=gene-tin

