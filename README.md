# tin
**Transcript Integrity Number (TIN)** calculation

TIN denotes how evenly a feature is covered by reads. By default, the script calculate the coverage evenness across a 
gene by considering the coverage across all exons.

In addition, the script also checks for *overlap transcription* where a transcript from the neighboring gene overlap with another gene.
It is important to take a note of such cases, as it can affect th differential expression results.

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

`--read_len` Read length. Default is 100 bp.

`--bg` If specified background noise will be subtracted
 
`--n` No. of bases to be subtracted from each ends of the feature to calculate effective length. Default is 50.
    
**How is tin calculated**

Overlapping features specified by `--feat` are merged together for a gene and tin is calculated 

TO Do: add formula  here.

**Transcript overlaps (run-ins)**

Transcript overlaps are flagged using `overlap.py` script.

	tincheck overlap --ann data/ann.gtf data/sample.bam >sample_overlap.txt

Genes are checked for overlap only if 

    gene-tin < tin-cutoff and gene-count > count-cutoff
    
 ie, genes with enough coverage but low tins are candidates for potential overlap.
 
If the count and tin values of neighboring genes and the inter-genic regions in the gene sense strand is 
comparable to the gene count and tin, then that gene is considered to have an overlap from neighboring transcripts.


More specifically, a gene is considered to have a neighboring transcript overlap, if **either of** the following conditions are satisfied.

    
    1. left gene tin in gene sense strand >= gene tin and left intergenic tin in gene sense strand >=gene-tin
    2. right gene tin in gene sense strand >= gene tin and right intergenic tin in gene sense strand >=gene-tin

