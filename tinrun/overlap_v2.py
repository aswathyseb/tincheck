import csv, sys, os, math, random
import argparse, subprocess, string

import plac
from tinrun import tin

"""
This script checks for run-in transcription overlap.
Inputs are bam files and gff3 file.

Dependencies:
    - featureCounts (conda install -c bioconda subread)
    - bedtools      (conda install -c bioconda bedtools)
    - samtools      (conda install -c bioconda samtools)
    
How it works?

Input(s) : bam file and annotation file
1. Intergenic and intronic annotation files (if background is specified) are created from the input ann file.
2. A bed file merging overlapping features are created (used to get depth per gene while calculating tin)
3. For each gene and intergenic region, counts are calculated in the sense and antisense strand.
4. Depending on gene strand, sense or antisense counts and tins are obtained and compared to that of the gene to
    determine the overlap.

"""

TMP = tin.TMP


class Measures:
    def __init__(self, count, tin):
        self.count = count
        self.tin = tin
        self.sense_count = ""
        self.sense_tin = ""
        self.antisense_count = ""
        self.antisense_tin = ""


def combine_counts(sense_counts, antisense_counts):
    """
    Sum up sense and antisense counts.
    Inputs are dictionary of list of dictionaries.
    Returns a dictionary of list of dictionaries.
    """
    total = dict()
    for key, vals in sense_counts.items():
        # print("***",key, vals)
        if key == "samples":
            total[key] = vals
            continue
        counts = []

        for idx, c in enumerate(vals):
            antisense_vals = antisense_counts[key]

            sample_sense_count = c
            sample_antisense_count = antisense_vals[idx]

            sample_count = sample_sense_count + sample_antisense_count

            counts.append(sample_count)

        total[key] = counts

    return total


def get_neighbors(genes, short):
    """
    Add left and right gene information to the dictionary.
    If a gene is completely inside another, left and right gene information is left empty.
    """

    left_gene, right_gene, = "", ""
    left_gene_strand, right_gene_strand = "", ""

    # Filter short genes
    filtered = {k: v for k, v in genes.items() if k not in short}

    vals = list(filtered.values())

    # Sort genes based on chrom and start
    vals.sort(key=lambda x: (x.chrom, x.start))

    for idx, gene in enumerate(vals):

        gid, chrom = gene.uid, gene.chrom

        # Not the very first or the very last genes.
        if idx != 0 and idx != len(vals) - 1:
            prev_gene = vals[idx - 1]
            next_gene = vals[idx + 1]
            prev_chrom, next_chrom = prev_gene.chrom, next_gene.chrom

            if chrom != next_chrom:
                right_gene, right_gene_strand = "", ""
                left_gene, left_gene_strand = prev_gene.uid, prev_gene.strand

            elif chrom != prev_chrom:
                left_gene, left_gene_strand = "", ""
                right_gene, right_gene_strand = next_gene.uid, next_gene.strand

            else:
                left_gene = prev_gene.uid
                left_gene_strand = prev_gene.strand
                right_gene = next_gene.uid
                right_gene_strand = next_gene.strand

        # first gene
        if idx == 0:
            next_gene = vals[idx + 1]
            left_gene, left_gene_strand = "", ""
            right_gene, right_gene_strand = next_gene.uid, next_gene.strand

        # last gene
        if idx == len(vals) - 1:
            right_gene, right_gene_strand = "", ""
            prev_gene = vals[idx - 1]
            left_gene, left_gene_strand = prev_gene.uid, prev_gene.strand

        gene.left_gene = left_gene
        gene.left_gene_strand = left_gene_strand
        gene.right_gene = right_gene
        gene.right_gene_strand = right_gene_strand
        genes[gid] = gene

    for k, v in genes.items():
        if k in short:
            v.left_gene, v.right_gene = "", ""
            v.left_gene_strand, v.right_gene_strand = "", ""
        genes[k] = v

    return genes


def bam_to_genome(bam):
    """
    Generate a genome file from bam header.
    The genome file has the format <chrom><tab>size>
    """
    refs = tin.parse_bam_header(bam)
    outfile = os.path.join(TMP, "chrom.txt")
    fh = open(outfile, "w")
    for r in refs:
        out = "\t".join([r[0], str(r[1])])
        fh.write(out + "\n")

    return outfile


def make_bed(coords, outfile="intergenic.bed"):
    """
    Create a bed file from a dictionary with name as key and coords as vals
    """
    outfile = os.path.join(TMP, outfile)
    fout = open(outfile, "w")

    for name, vals in coords.items():
        chrom, start, end = vals[0], vals[1], vals[2]
        out = "\t".join([chrom, str(start), str(end), name, ".", "+"])
        tin.write_to_file([out], fout)
    return outfile


def make_intergenic_files(coords):
    ig_ann = make_intergenic_gtf(coords)
    ig_bed = make_bed(coords)
    return ig_bed, ig_ann


def make_intergenic_gtf(coords, outfile="intergenic.gtf"):
    """
    Create a minimal gtf file for intergenic regions.
    """

    outfile = os.path.join(TMP, outfile)
    fout = open(outfile, "w")

    idx = 0
    for key, vals in coords.items():
        gid = key
        chrom, start, end = vals[0], vals[1], vals[2]
        feat_ln = tin.make_gtf_feature_rows(chrom, start, end, "+", gid, idx, "intergenic")
        tin.write_to_file([feat_ln], fout)
    return outfile


def name_intergenic_regions(intergenic, genes, short):
    """"
    Name the intergenic intervals with the geneids of the neighboring genes in the format gene_id1:geneid2
    if (igstart, igend) are in between gene1-start and gene2-start a, then it will get the name gene1:gene2.
    Returns a dictionary with intergenic-name as key and its coordinates as values.
    """
    ignames = dict()

    gcoords, igcoords, gnames = [], [], []

    # Filter short genes
    filtered = {k: v for k, v in genes.items() if k not in short}

    # Get gene (chrom, start,end, gid) in list
    for gene, vals in filtered.items():
        chrom, start, end, gid = vals.chrom, vals.start, vals.end, vals.uid
        gcoords.append((chrom, int(start), int(end), gid))
        gnames.append(gid)

    # Get the intergenic (chrom,start,end in a list)
    for item in intergenic:
        chrom, start, end = item.split("\t")
        start, end = int(start) + 1, int(end)
        igcoords.append((chrom, start, end))

    # Sort gene list and intergenic list by chrom and start
    gcoords.sort(key=lambda x: (x[0], x[1]))
    igcoords.sort(key=lambda x: (x[0], x[1]))

    # Go through the intergenic list and check if its coordinates are between the neighboring gene starts.
    for i in igcoords:
        chrom, start, end = i

        for idx, g in enumerate(gcoords):

            if idx == len(gcoords) - 1:
                # Last one
                gchrom1, gstart1, gend1, gid1 = gcoords[idx - 1]
                gchrom2, gstart2, gend2, gid2 = g
                if chrom == gchrom1 and chrom == gchrom2 and start > gstart1 and end < gstart2:
                    name = ":".join([gid1, gid2])
                    ignames[name] = i
                break

            gchrom1, gstart1, gend1, gid1 = g
            gchrom2, gstart2, gend2, gid2 = gcoords[idx + 1]
            if chrom == gchrom1 and chrom == gchrom2 and start > gstart1 and end < gstart2:
                name = ":".join([gid1, gid2])
                ignames[name] = i
                break
    return ignames


def get_intergenic_coords(genes, genome):
    """
    Create intergenic regions using bedtools.
    """

    genes_bed = os.path.join(TMP, "genes.bed")

    fh = open(genes_bed, "w")

    # Get all gene coordinates.
    for k, v in genes.items():
        gene = genes[k]
        chrom, start, end = gene.chrom, gene.start - 1, gene.end
        fh.write("\t".join([chrom, str(start), str(end)]))
        fh.write("\n")

    fh.close()

    try:
        cmd = f"cat {genes_bed}| bedtools sort -faidx {genome} -i - | bedtools complement -i - -g {genome}"
        res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, check=True, universal_newlines=True)
        # print(cmd)

    except:
        print("Error creating the intergenic intervals.")
        sys.exit()

    res_list = res.stdout.split("\n")
    res_list = res_list[:-1]

    return res_list


def get_short_genes(genes):
    """
    Filter short genes that are completely inside another.
    """
    short = []

    # get short genes
    for gene1, val1 in genes.items():

        chrom1, start1, end1, gid1 = val1.chrom, val1.start, val1.end, val1.uid

        for gene2, val2 in genes.items():
            chrom2, start2, end2, gid2 = val2.chrom, val2.start, val2.end, val2.uid

            if chrom1 != chrom2:
                continue

            if start1 >= start2 and end1 <= end2 and gid1 != gid2:
                short.append(gid1)

    return short


def get_sample_vals(info, sample_index):
    """
    Extract sample specific count from a master dictionary of lists.
    """
    sample_vals = dict()
    for key, vals in info.items():
        if key == "samples":
            continue
        sample_vals[key] = vals[sample_index]

    return sample_vals


def collect_measures(counts, tins, sense_counts={}, antisense_counts={}, sense_tins={},
                     antisense_tins={}, strand="both"):
    """
    Collect exp,tin, obs_tin, tin_sense, tin_antisense in a single dictionary.
    """
    gtins = dict()

    for gene, count in counts.items():
        tin = tins[gene][2]
        m = Measures(count=count, tin=tin)

        if strand != "both":
            sense_count = sense_counts[gene]
            antisense_count = antisense_counts[gene]
            sense_tin = sense_tins[gene][2]
            antisense_tin = antisense_tins[gene][2]

            m.sense_tin = sense_tin
            m.sense_count = sense_count
            m.antisense_tin = antisense_tin
            m.antisense_count = antisense_count

        gtins[gene] = m

    return gtins


def get_neighbor_values(lgene, rgene, ileft, iright, gmeasures, igmeasures):
    """
    Collect the tin and count values of left gene, left  intergenic region, right gene and right intergenic region
    in unstranded mode.
    """

    ileft_tin = igmeasures[ileft].tin if ileft in igmeasures else None
    iright_tin = igmeasures[iright].tin if iright in igmeasures else None
    ileft_count = igmeasures[ileft].count if ileft in igmeasures else None
    iright_count = igmeasures[iright].count if iright in igmeasures else None

    gleft_tin = gmeasures[lgene].tin
    gright_tin = gmeasures[rgene].tin
    gleft_count = gmeasures[lgene].count
    gright_count = gmeasures[rgene].count

    ig_vals = (ileft_tin, iright_tin, ileft_count, iright_count)
    gene_vals = (gleft_tin, gright_tin, gleft_count, gright_count)

    return gene_vals, ig_vals


def check_runin(data, gmeasures, igmeasures, strand, tin_cutoff=40, count_cutoff=40):
    """
    Check runins (overlap) for genes with tin <= tin_cutoff and count >=count_cutoff.
    Tin  and count values of the feature is checked with the tin and count values in the feature strand
    of the neighboring genes on both sides.
    """

    ileft_tin, gleft_tin, iright_tin, gright_tin = 0, 0, 0, 0
    ileft_count, gleft_count, iright_count, gright_count = 0, 0, 0, 0

    runins, params = dict(), dict()

    for gene, vals in gmeasures.items():

        gene_tin = vals.tin
        gene_count = vals.count

        # Checking for overlap only for genes with tin in range 1 to 40 and (count!=0 and count >count_cutoff).
        # ie, enough reads/expression but not good tin.
        if gene_tin > tin_cutoff or gene_tin == 0 or gene_count < count_cutoff:
            continue

        gene_strand = data[gene].strand
        lgene_strand = data[gene].left_gene_strand
        rgene_strand = data[gene].right_gene_strand
        lgene = data[gene].left_gene
        rgene = data[gene].right_gene

        # The very first, the very last and short genes.
        if lgene == "" or rgene == "":
            continue

        # Intergenic dict keys : left_right
        ileft = ":".join([lgene, gene])
        iright = ":".join([gene, rgene])

        # Unstranded mode check.
        if strand == "both":
            gene_vals, ig_vals = get_neighbor_values(lgene, rgene, ileft, iright, gmeasures, igmeasures)
            gleft_tin, gright_tin, gleft_count, gright_count = gene_vals
            ileft_tin, iright_tin, ileft_count, iright_count = ig_vals

        # Stranded mode check.
        else:

            if lgene_strand == "+":
                gleft_tin, gleft_count = gmeasures[lgene].sense_tin, gmeasures[lgene].sense_count
                ileft_tin = igmeasures[ileft].sense_tin if ileft in igmeasures else None
                ileft_count = igmeasures[ileft].sense_count if ileft in igmeasures else None

            if lgene_strand == "-":
                gleft_tin, gleft_count = gmeasures[lgene].antisense_tin, gmeasures[lgene].antisense_count
                ileft_tin = igmeasures[ileft].sense_tin if ileft in igmeasures else None
                ileft_count = igmeasures[ileft].sense_count if ileft in igmeasures else None

            if  rgene_strand == "+":
                gright_tin, gright_count = gmeasures[rgene].antisense_tin, gmeasures[rgene].antisense_count
                iright_tin = igmeasures[iright].antisense_tin if iright in igmeasures else None
                iright_count = igmeasures[iright].antisense_count if iright in igmeasures else None

            if rgene_strand == "-":
                gright_tin, gright_count = gmeasures[rgene].sense_tin, gmeasures[rgene].sense_count
                iright_tin = igmeasures[iright].antisense_tin if iright in igmeasures else None
                iright_count = igmeasures[iright].antisense_count if iright in igmeasures else None



                    # if gene_strand == "+" and strand == "antisense":
            #     ileft_tin = igmeasures[ileft].antisense_tin if ileft in igmeasures else None
            #     iright_tin = igmeasures[iright].antisense_tin if iright in igmeasures else None
            #     ileft_count = igmeasures[ileft].antisense_count if ileft in igmeasures else None
            #     iright_count = igmeasures[iright].antisense_count if iright in igmeasures else None
            #
            #     gleft_tin = gmeasures[lgene].antisense_tin if lgene_strand == "+" else gmeasures[lgene].sense_tin
            #     gright_tin = gmeasures[rgene].antisense_tin if rgene_strand == "+" else gmeasures[rgene].sense_tin
            #     gleft_count = gmeasures[lgene].antisense_count if lgene_strand == "+" else gmeasures[lgene].sense_count
            #     gright_count = gmeasures[rgene].antisense_count if rgene_strand == "+" else gmeasures[rgene].sense_count
            #
            # if gene_strand == "-" and strand == "antisense":
            #     ileft_tin = igmeasures[ileft].sense_tin if ileft in igmeasures else None
            #     iright_tin = igmeasures[iright].sense_tin if iright in igmeasures else None
            #     ileft_count = igmeasures[ileft].sense_count if ileft in igmeasures else None
            #     iright_count = igmeasures[iright].sense_count if iright in igmeasures else None
            #
            #     gleft_tin = gmeasures[lgene].antisense_tin if lgene_strand == "-" else gmeasures[lgene].sense_tin
            #     gright_tin = gmeasures[rgene].antisense_tin if rgene_strand == "-" else gmeasures[rgene].sense_tin
            #     gleft_count = gmeasures[lgene].antisense_count if lgene_strand == "-" else gmeasures[lgene].sense_count
            #     gright_count = gmeasures[rgene].antisense_count if rgene_strand == "-" else gmeasures[rgene].sense_count
            #
            # if gene_strand == "+" and strand == "sense":
            #     ileft_tin = igmeasures[ileft].sense_tin if ileft in igmeasures else None
            #     iright_tin = igmeasures[iright].sense_tin if iright in igmeasures else None
            #     ileft_count = igmeasures[ileft].sense_count if ileft in igmeasures else None
            #     iright_count = igmeasures[iright].sense_count if iright in igmeasures else None
            #
            #     gleft_tin = gmeasures[lgene].sense_tin if lgene_strand == "+" else gmeasures[lgene].antisense_tin
            #     gright_tin = gmeasures[rgene].sense_tin if rgene_strand == "+" else gmeasures[rgene].antisense_tin
            #     gleft_count = gmeasures[lgene].sense_count if lgene_strand == "+" else gmeasures[lgene].antisense_count
            #     gright_count = gmeasures[rgene].sense_count if rgene_strand == "+" else gmeasures[rgene].antisense_count
            #
            # if gene_strand == "-" and strand == "sense":
            #     ileft_tin = igmeasures[ileft].antisense_tin if ileft in igmeasures else None
            #     iright_tin = igmeasures[iright].antisense_tin if iright in igmeasures else None
            #     ileft_count = igmeasures[ileft].antisense_count if ileft in igmeasures else None
            #     iright_count = igmeasures[iright].antisense_count if iright in igmeasures else None
            #
            #     gleft_tin = gmeasures[lgene].sense_tin if lgene_strand == "-" else gmeasures[lgene].antisense_tin
            #     gright_tin = gmeasures[rgene].sense_tin if rgene_strand == "-" else gmeasures[rgene].antisense_tin
            #     gleft_count = gmeasures[lgene].sense_count if lgene_strand == "-" else gmeasures[lgene].antisense_count
            #     gright_count = gmeasures[rgene].sense_count if rgene_strand == "-" else gmeasures[rgene].antisense_count

        # If the intergenic region is None, set the values to be equivalent to gene,
        # so that a general condition can be applied.

        ileft_tin = gene_tin if ileft_tin is None else ileft_tin
        iright_tin = gene_tin if iright_tin is None else iright_tin
        ileft_count = gene_count if ileft_count is None else ileft_count
        iright_count = gene_count if iright_count is None else iright_count

        if left_tin_condition(gene_tin, gleft_tin, ileft_tin) and \
                left_count_condition(gene_count, gleft_count, ileft_count) and \
                right_tin_condition(gene_tin, gright_tin, iright_tin) and \
                right_count_condition(gene_count, gright_count, iright_count):
            runins[gene] = f"both"
        elif left_tin_condition(gene_tin, gleft_tin, ileft_tin) and \
                left_count_condition(gene_count, gleft_count, ileft_count):
            runins[gene] = "left"
        elif right_tin_condition(gene_tin, gright_tin, iright_tin) and \
                right_count_condition(gene_count, gright_count, iright_count):
            runins[gene] = "right"

    return runins


def left_tin_condition(gene_tin, gleft_tin, igleft_tin):
    return gleft_tin >= gene_tin and igleft_tin >= gene_tin


def right_tin_condition(gene_tin, gright_tin, igright_tin):
    return gright_tin >= gene_tin and igright_tin >= gene_tin


def left_count_condition(gene_count, gleft_count, igleft_count):
    return gleft_count >= gene_count and igleft_count >= gene_count


def right_count_condition(gene_count, gright_count, igright_count):
    return gright_count >= gene_count and igright_count >= gene_count


@plac.pos('bams', "comma separated bam files.")
@plac.opt('ann', help="annotation file in GTF or GFF3 format.")
@plac.opt('feat', type=str, help="feature in the 3rd column of the anntation file on which TIN needs to be calculated")
@plac.opt('groupby', type=str, help="attribute by which features need to be combined , eg: gene_id")
@plac.opt('strand', type=str, help="strand on which tin should be calculated, eg: both, sense or antisense")
@plac.opt('libtype', type=str, help="library type eg: paired or single")
@plac.flg('bg', help="when specified background noise will be subtracted")
@plac.opt('n', type=int, help="""number of bases to be subtracted from
                              each ends of the feature to calculate effective length""")
@plac.opt('tin_cutoff', type=int, help="""tin cutoff to be used for checking overlaps 
                    Genes with tin < tin-cutoff and count >count_cutoff are checked for overlaps""")
@plac.opt('count_cutoff', type=int,
          help="count cutoff to be used for checking overlaps. ie, minimum reads required to be considered a gene as being expressed")
def run(bams, ann="", feat='exon', groupby='gene_id', strand='both', libtype='single', bg=False, n=50,
        tin_cutoff=40, count_cutoff=40):
    bg_file, intron_len, intron_gtf = None, None, None

    # Check if inputs are valid.
    tin.check_inputs(bams, strand, libtype, ann)

    # Create TMP directory
    os.makedirs(TMP)

    # Create run-log.
    log = os.path.join(TMP, "runlog.txt")

    # Make a string of bam files.
    bams = bams.split(',')
    bam = " ".join(bams)

    # Extract features from annotation file.
    genes, features = tin.collect_features(ann=ann, feat_type=feat, groupby=groupby)

    # Get genes that are completely inside another.
    short = get_short_genes(genes)

    # Add neighbors
    genes = get_neighbors(genes=genes, short=short)

    # Merge overlapping features by groupby attribute.
    merged = tin.merge_features(features)

    # Make a bed file of merged features.
    bed = tin.make_feature_bed(genes, merged)

    # Generate genome file.
    genomefile = bam_to_genome(bam)

    # Get intergenic regions
    interg_coords = get_intergenic_coords(genes=genes, genome=genomefile)

    # Name intergenic regions
    intergenic = name_intergenic_regions(interg_coords, genes, short)

    # Make intergenic bed and gff files.
    ig_bed, ig_ann = make_intergenic_files(intergenic)

    # Get intron coordinates for background noise calculation.
    if bg:
        introns = tin.get_intron_coords(merged)
        intron_len = tin.get_effective_length(introns, size=0)
        intron_gtf = tin.make_intron_gtf(genes, introns)

    # Paired or single.
    paired = True if libtype == "paired" else False

    # Get gene counts in both sense and antisense.
    gcounts_sense = tin.get_counts(bam=bam, ann=ann, strand="sense", feat=feat,
                                   groupby=groupby, paired=paired, flag="gene_sense", log=log)

    gcounts_antisense = tin.get_counts(bam=bam, ann=ann, strand="antisense", feat=feat,
                                       groupby=groupby, paired=paired, flag="gene_antisense", log=log)

    # Get intergenic counts in both sense and antisense.
    igcounts_sense = tin.get_counts(bam=bam, ann=ig_ann, strand="sense", feat="intergenic",
                                    groupby=groupby, paired=paired, flag="ig_sense", log=log)

    igcounts_antisense = tin.get_counts(bam=bam, ann=ig_ann, strand="antisense", feat="intergenic",
                                        groupby=groupby, paired=paired, flag="ig_antisense", log=log)

    # Get the observed gene counts and intergenic according to the library strand.
    if strand == "both":
        # Gene counts in unstranded mode. Sum up both sense and antisense counts.
        gene_counts = combine_counts(gcounts_sense, gcounts_antisense)
        ig_counts = combine_counts(igcounts_sense, igcounts_antisense)
    else:
        gene_counts = gcounts_antisense if strand == "antisense" else gcounts_sense
        ig_counts = igcounts_antisense if strand == "antisense" else igcounts_sense

    # Get feature length
    feat_len = tin.get_effective_length(merged, n)

    # Calculate expected tin.
    # exp_tins = tin.get_exp_tin(counts=gene_counts, paired=paired, read_len=read_len, feat_len=feat_len)

    # Sample specific TIN calculations begins.
    runins, obs_tins, exp_tins = dict(), dict(), dict()

    for idx, bam in enumerate(bams):

        # Parse out the gene counts and intergenic counts for the current bam.
        sample_gene_counts = get_sample_vals(info=gene_counts, sample_index=idx)
        sample_ig_counts = get_sample_vals(info=ig_counts, sample_index=idx)

        sample_gcounts_sense = get_sample_vals(info=gcounts_sense, sample_index=idx)
        sample_gcounts_antisense = get_sample_vals(info=gcounts_antisense, sample_index=idx)

        sample_ig_counts_sense = get_sample_vals(info=igcounts_sense, sample_index=idx)
        sample_ig_counts_antisense = get_sample_vals(info=igcounts_antisense, sample_index=idx)

        # Calculate background noise if "--bg" is specified.
        if bg:
            bg_file = tin.get_background_noise(bam=bam, strand=strand,
                                               intron_ann=intron_gtf, intron_len=intron_len, groupby=groupby,
                                               paired=paired)

        # Extract primary alignments
        pbam = tin.get_primary_aln(bam)

        if strand == "both":
            # Get gene tin.
            gene_tin = tin.get_gene_tin(bam=pbam, bed=bed, strand=strand, bgfile=bg_file, size=n, flag="gene")

            # Get intergenic tin.
            ig_tin = tin.get_gene_tin(bam=pbam, bed=ig_bed, strand=strand, size=0, flag="ig")

            # Get the gene counts and tins calculated for the sample in a single place.
            gene_measures = collect_measures(counts=sample_gene_counts, tins=gene_tin,
                                             strand=strand)

            # Get intergenic counts and tins for the sample single place.
            intergenic_measures = collect_measures(counts=sample_ig_counts, tins=ig_tin, strand=strand)

        else:

            # Get gene sense and antisense tin.
            gtin_sense = tin.get_gene_tin(bam=pbam, bed=bed, strand="sense", bgfile=bg_file, size=n, flag="gene_sense")
            gtin_antisense = tin.get_gene_tin(bam=pbam, bed=bed, strand="antisense", bgfile=bg_file, size=n,
                                              flag="gene_antisense")

            # Get intergenic sense and antisense tin. Intergenic bed has + as sense strand.
            igtin_sense = tin.get_gene_tin(bam=pbam, bed=ig_bed, strand="sense", size=0, flag="ig_sense")
            igtin_antisense = tin.get_gene_tin(bam=pbam, bed=ig_bed, strand="antisense", size=0, flag="ig_antisense")

            # Get tin according to the library strand.
            gene_tin = gtin_antisense if strand == "antisense" else gtin_sense
            ig_tin = igtin_antisense if strand == "sense" else igtin_sense

            # Get the counts and the different types of tins calculated for a gene in a single place.
            gene_measures = collect_measures(counts=sample_gene_counts, tins=gene_tin,
                                             sense_counts=sample_gcounts_sense,
                                             antisense_counts=sample_gcounts_antisense,
                                             sense_tins=gtin_sense, antisense_tins=gtin_antisense, strand=strand)

            # Get intergenic counts and tins in a single place.
            intergenic_measures = collect_measures(counts=sample_ig_counts, tins=ig_tin,
                                                   sense_counts=sample_ig_counts_sense,
                                                   antisense_counts=sample_ig_counts_antisense,
                                                   sense_tins=igtin_sense,
                                                   antisense_tins=igtin_antisense, strand=strand)

        # Check overlap
        overlaps = check_runin(data=genes, gmeasures=gene_measures, igmeasures=intergenic_measures,
                               strand=strand, tin_cutoff=tin_cutoff, count_cutoff=count_cutoff)

        # Collect obs_tins.
        for uid, vals in gene_tin.items():
            if uid == "samples":
                obs_tins.setdefault(uid, []).append(vals[0])
                exp_tins.setdefault(uid, []).append(vals[1])
            else:
                obs_tins.setdefault(uid, []).append(vals[2])
                exp_tins.setdefault(uid, []).append(vals[3])

        # Collect runins
        for gene in gene_counts:
            if gene == "samples":
                continue
            overlap = overlaps[gene] if gene in overlaps else ""
            runins.setdefault(gene, []).append(overlap)
        runins.setdefault('samples', []).append(tin.get_filename(bam) + "_overlap")

    # Collect all results, gid,counts, exp_tin, obs_tin
    results = dict()

    for k, v in gene_counts.items():
        vals = list()
        vals.append(feat_len[k])
        vals.extend(v)
        vals.extend(exp_tins[k])
        vals.extend(obs_tins[k])
        vals.extend(runins[k])
        results[k] = vals

    # Print results
    for gene, vals in results.items():
        gene = "target_id" if gene == "samples" else gene
        vals = map(str, vals)
        vals = "\t".join(vals)
        out = "\t".join([gene, vals])
        print(out)

    # Clean up temporary files.
    cmd = f'rm -rf {TMP}'
    # os.system(cmd)


if __name__ == "__main__":
    run()
