import csv, sys, os, math, random, re
import argparse, subprocess, string
import uuid
import plac

"""
This script produces raw read-count counts, expected-tins and observed-tins.
Inputs are bam files and/or gtf file.

Dependencies:
    - featureCounts (conda install -c bioconda subread)
    - bedtools      (conda install -c bioconda bedtools)
    - samtools      (conda install -c bioconda samtools)
    
Program calculates tin values at a specified level (eg: exon, CDS etc grouped by gene or transcript )
It can also calculate TINs from bam file when there is no annotation file.
This is useful for calculating transcript tins from transcriptome bam file.

"""

TMP = "tmp_" + str(uuid.uuid4())[:3]


class Features:
    def __init__(self, chrom, start, end, strand, uid, type):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.type = type
        self.uid = uid
        self.left_gene = ''
        self.right_gene = ''
        self.left_gene_strand = ''
        self.right_gene_strand = ''

    def get_params(self):
        return self.chrom, self.start, self.end, self.strand, self.uid, self.type


def get_filename(x):
    return os.path.splitext(os.path.basename(x))[0]


def add_string(elms, string):
    def attach(x):
        return x + string

    elms = list(map(attach, elms))
    return elms


def store_counts(fname):
    """
    Store featureCounts output in a dictionary.
    """
    store = dict()
    stream = csv.reader(open(fname, 'rt'), delimiter="\t")

    # skip header rows.
    next(stream)

    header = next(stream)
    samples = header[6:]
    samples = map(get_filename, samples)
    samples = add_string(samples, "_count")
    store.setdefault('samples', []).extend(samples)

    for row in stream:
        # Skip first 5 columns as counts start at 6th column
        gid, glen, cols = row[0], row[5], row[6:]
        for count in cols:
            store.setdefault(gid, []).append(int(count))

    return store


def get_counts(bam, ann, strand="", feat="exon", groupby="gene", paired=False, flag="", log="runlog.txt", ):
    """
    Use featureCounts to get the raw counts.
    """

    countsfile = f"counts_{flag}.txt" if len(bam.split(' ')) > 1 else get_filename(bam) + "_counts.txt"
    countsfile = countsfile.replace("_counts.txt", f'_{flag}_counts.txt') if flag else countsfile
    countsfile = os.path.join(TMP, countsfile)

    # Set paired end or single end.
    pflag = '-p' if paired else ""

    group_attr = "gene_id" if feat == "intergenic" else groupby

    # Set strandedness.
    if strand == "both":
        sflag = ""
    else:
        sflag = '-s 2' if strand == "sense" else '-s 1'

    # Make the command.
    cmd = f'featureCounts --primary -M {pflag} {sflag} -t {feat} -O -g {group_attr} -a {ann} -o {countsfile} {bam} 2>>{log}'

    # print(cmd)

    # Run command and get counts.
    out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
                         universal_newlines=True)

    # Exit on error.
    exit_on_error(out.stderr) if out.stderr and out.returncode != 0 else None

    # Store counts in a dictionary.
    counts = store_counts(countsfile)

    return counts


def get_intron_coords(store):
    """
    Collect all intron coordinates.
    Input store is a dictionary of list where each element in the list is (start,end) of exons/cds
    Returns a dictionary of list where each element in the list is (start,end) of intron.
    """
    introns = dict()
    for key, vals in store.items():

        if len(vals) == 1:  # single exon genes
            continue
        # Collect intron starts and ends.
        istarts = [end + 1 for (start, end) in vals[:-1]]
        iends = [start - 1 for (start, end) in vals[1:]]
        coords = list(zip(istarts, iends))
        introns.setdefault(key, []).extend(coords)
    return introns


def make_intron_gtf(genes, introns):
    """
    Create a GTF file for introns.
    """
    # intron gtf file.
    intron_ann = os.path.join(TMP, "introns.gtf")
    fi = open(intron_ann, "w")

    for gene, coords in introns.items():
        # Get gene coords
        gene_vals = genes[gene]
        chrom, start, end, strand = gene_vals.chrom, gene_vals.start, gene_vals.end, gene_vals.strand
        gid, gtype = gene_vals.uid, gene_vals.type

        # Write gene line.
        gene_attr = f'ID \"{gid}\";gene_id \"{gid}\"'
        parent_ln = make_ann_row(chrom, gtype, start, end, strand, gene_attr)
        fi.write(parent_ln + "\n")

        for idx, coord in enumerate(coords):
            # Write intron line.
            istart, iend = str(coord[0]), str(coord[1])
            intron_ln = make_gtf_feature_rows(chrom, istart, iend, strand, gid, idx=idx, feat="intron")
            fi.write(intron_ln + "\n")

    return intron_ann


def get_background_noise(bam, strand, intron_ann, intron_len, groupby, paired):
    """
    Calculate background noise from intron specific gff file.
    Returns a dictionary with background noise.
    """
    noise_file = get_filename(bam) + "_noise.txt"
    noise_file = os.path.join(TMP, noise_file)
    nout = open(noise_file, "w")

    intron_counts = get_counts(bam=bam, ann=intron_ann, strand=strand,
                               feat="intron", groupby=groupby, paired=paired, flag="intron")

    for key, vals in intron_counts.items():
        if key == "samples":
            continue
        for val in vals:
            c, l = val, intron_len.get(key, 0)
            noise = 0 if c == 0 or l == 0 else round(int(c) / int(l), 4)
            nout.write(f'{key}\t{str(noise)}')
            nout.write("\n")

    return noise_file


def extract_attr(attr_str, attr_tag):
    attrs = attr_str.split(";")
    attrs = [a.strip() for a in attrs]  # remove all spaces
    try:
        feat_id = [val for val in attrs if val.startswith(attr_tag)][0]
    except:
        print(f'{attr_tag} not found for the feature in annotation file')
        sys.exit()
    feat_id = feat_id.replace(attr_tag + "=", "")  # for gff
    feat_id = feat_id.replace(attr_tag + " ", "")  # for gtf
    feat_id = feat_id.replace('"', "")  # remove all quotes.
    # feat_id = feat_id.split(",")  # for cases like Parent=mRNA1,mRNA2,mRNA3 in gff3
    return feat_id


def parse_row(row):
    # chrom, start, end, strand, feat
    return row[0], int(row[3]), int(row[4]), row[6], row[2], row[8]


def merge(elms):
    overlapped = list()
    elms = sorted(elms)
    elms.append((1000000000, 1000000000))

    start = end = 0

    for (x1, y1), (x2, y2) in zip(elms, elms[1:]):

        start = x1 if start == 0 else start
        end = y1 if end == 0 else end

        # Case1 : left overlap
        if x2 <= end and y2 > end:
            end = y2

        # Case2 : complete overlap
        if x2 <= end and y2 <= end:
            start = start
            end = end

        # Case3 : no overlap
        if x2 > end:
            # print(start, end)
            overlapped.append((start, end))
            start = end = 0
    return overlapped


def merge_features(feats):
    """
    Make a union of features for each gene by merging overlapping features.
    """
    merged = dict()
    for gene, elms in feats.items():
        combined = merge(elms)
        merged.setdefault(gene, combined)
    return merged


def collect_features(ann, feat_type="exon", groupby="gene"):
    """
    Collect the features (exons or cds)  grouped by 'groupby'  attribute in a dictionary.
    """
    feat_store, genes = dict(), dict()

    # Parse file and collect features.
    stream = csv.reader(open(ann), delimiter="\t")

    for row in stream:
        # Continue if header.
        if row[0].startswith("#"):
            continue

        chrom, start, end, strand, feat, attrs = parse_row(row)

        # Make sure start < end
        if start > end:
            start, end = end, start

        if feat == feat_type:
            gid = extract_attr(attrs, groupby)

            genes[gid] = Features(chrom=chrom, start=0, end=0, strand=strand, uid=gid, type='')
            feat_store.setdefault(gid, []).append((start, end))

    # Get gene start and end from the collected exons.
    for k, v in feat_store.items():
        first, last = sorted(v)[0], sorted(v)[-1]
        start, end = first[0], last[1]
        genes[k].start = start
        genes[k].end = end

    # Exit if feat_type or group_by not in the annotation file.
    empty_feat_store = [k for k, a in feat_store.items() if a == []]

    if not genes or len(empty_feat_store) == len(genes):
        print(f'{feat_type} do not match with 3rd column in annotation file. Exiting.')
        sys.exit()
    return genes, feat_store


def get_ref_names(string):
    """
    Returns the reference sequences names and length from the sam header string.
    """
    if string.startswith("@SQ"):
        elms = string.split("\t")
        elms = [e.strip() for e in elms]
        name = [e for e in elms if e.startswith("SN:")][0]
        slen = [e for e in elms if e.startswith("LN:")][0]
        name = name.replace("SN:", "")
        slen = slen.replace("LN:", "")

        return name, int(slen)


def parse_bam_header(bam):
    """
    Parses bam header to get reference names and lengths.
    """
    # Get one bam file if there are more than one.
    vals = bam.split(' ')
    aln = vals[0] if len(vals) > 1 else bam

    # Extract header from bam file.
    cmd = f'samtools view -H {aln}'
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, check=True, universal_newlines=True)

    # Parse header to extract sequence name and length.
    res_list = res.stdout.split("\n")
    refs = map(get_ref_names, res_list)
    refs = filter(None, refs)
    return refs


def bam_to_gtf(bam):
    refs = parse_bam_header(bam)
    start, strand = 1, "+"
    outfile = os.path.join(TMP, "ann.gtf")
    fh = open(outfile, "w")
    for r in refs:
        chrom, end = r[0], r[1]
        gene_ln = make_gtf_feature_rows(chrom, start, end, strand, chrom, idx=0, feat="gene")
        trans_ln = make_gtf_feature_rows(chrom, start, end, strand, chrom, idx=0, feat="transcript")
        exon_ln = make_gtf_feature_rows(chrom, start, end, strand, chrom, idx=0, feat="exon")
        write_to_file([gene_ln, trans_ln, exon_ln], fh)

    return outfile


def make_ann_row(chrom, feat, start, end, strand, attr):
    return "\t".join([chrom, ".", feat, str(start), str(end), ".", strand, ".", attr])


def make_gtf_feature_rows(chrom, start, end, strand, gene_id, idx=0, feat="exon"):
    """
    Creates gtf feature rows from chr,start,end,strand,gene_id
    """
    # Make ids for gene and feature.
    num = str(idx + 1)
    # trans_id = gene_id + "_t" + num
    feat_id = gene_id + "_" + feat + num

    # Make 9th column
    attr = f'ID \"{feat_id}\";gene_id \"{gene_id}\";transcript_id \"{gene_id}\"'

    # Create feature rows.
    feat_ln = make_ann_row(chrom, feat, start, end, strand, attr)
    return feat_ln


def write_to_file(rows, fh):
    for row in rows:
        fh.write(row + "\n")


def make_feature_bed(genes, merged):
    """
    Makes a bed file of features.
    """
    bedfile = os.path.join(TMP, "features.bed")
    fh = open(bedfile, "w")

    for key, vals in genes.items():
        chrom, start, end, strand, gtype, gid = vals.chrom, vals.start, vals.end, vals.strand, vals.type, vals.uid

        # If there are no features for a gene, continue
        if not merged[gid]:
            continue

        coords = merged[gid]

        for c in coords:
            start, end = str(c[0] - 1), str(c[1])
            out = "\t".join([chrom, start, end, gid, ".", strand])
            fh.write(out + "\n")

    return bedfile


def get_depth(bam, bed, strand='sense', flag=""):
    """
    Use bedtools to get the depth at every position.
    """

    def get_strand_flag(strand):
        if strand == "both":
            lib_strand = ""
        else:
            lib_strand = '-S' if strand == "sense" else '-s'
        return lib_strand

    # Set strandedness.
    lib_strand = get_strand_flag(strand)

    # Output file.
    suffix1 = "depth.txt"
    outfile = get_filename(bam)
    suffix = f"_{flag}_{suffix1}" if flag else f"_{suffix1}"
    outfile = outfile.replace("_primary", "") + suffix
    outfile = os.path.join(TMP, outfile)

    # Command to get the depth.
    cmd2 = f'bedtools coverage -d -a {bed} -b {bam} {lib_strand} | cut -f 4,8 > {outfile}'
    # print(cmd2)

    d = subprocess.run(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
                       universal_newlines=True)
    # Exit on error.
    exit_on_error(d.stderr) if d.stderr and d.returncode != 0 else None

    return outfile


def get_gene_tin(bam, bed, strand, bgfile="", size=50, flag=""):
    # Get depth for each position in each gene.
    depth_file = get_depth(bam, bed, strand, flag=flag)

    # Get tin.
    tin_store = get_tin(depth_file, bgfile, size)

    # Collect sample name.
    sample = get_filename(bam).replace("_primary", "")
    # tin_store['samples'] = sample + "_obs_tin"
    sample_cols = [sample + "_obs_tin", sample + "_exp_tin"]
    tin_store.setdefault('samples', []).extend(sample_cols)
    return tin_store


def get_tin(depth_file, bgfile="", size=50):
    """
    Calculate tin from the depth file.
    Return effective_length, background_noise and tin_score.
    """

    covs, background, store = dict(), dict(), dict()

    curr, uid, path = "", "", ""

    if bgfile:
        background = store_background(bgfile)

    stream = csv.reader(open(depth_file), delimiter="\t")

    for row in stream:

        uid, depth = row[0], float(row[1])

        # The very first time
        if not curr:
            curr = uid

        if uid != curr:
            fid, tlen, bg, obs_tin, exp_tin = calculate_tin(coverages=covs, background=background, uid=curr, size=size)
            store[fid] = (tlen, bg, obs_tin, exp_tin)
            curr = uid
            covs = dict()

        covs.setdefault(uid, []).append(depth)

    # Print the last element
    fid, tlen, bg, obs_tin, exp_tin = calculate_tin(coverages=covs, background=background, uid=uid, size=size)

    store[fid] = (tlen, bg, obs_tin, exp_tin)

    return store


def calculate_tin(coverages, background, uid, size):
    """
    Calculates tin for the coverages and store in a dictionary
    where values are a tuple with(effective_len, background,tin-score)
    """

    vals = coverages[uid]

    # Get background.
    bg = float(background.get(uid, 0.0))

    # Subtract background.
    vals = subtract_background(vals, bg)

    # Make coverages float.
    vals = list(map(float, vals))

    # Calculate tin score.
    obs_tin, exp_tin, tlen = tin_score(cvg=vals, size=size)

    obs_tin = round(obs_tin, 1)
    exp_tin = round(exp_tin, 1)

    return uid, tlen, bg, obs_tin, exp_tin


def shannon_entropy(vals):
    """
    Calculate shannon's H = -sum(P*log(P)). Argument is a list of float numbers.
    """

    val_sum = sum(vals)
    entropy = 0.0

    for i in vals:
        entropy += (i / val_sum) * math.log(i / val_sum)

    return 0 if entropy == 0.0 else -entropy


def tin_score(cvg, size=50):
    """
    Calculate TIN score.
    cvg : coverage at each base position as a list.
    size : no. of bases to be omitted from the beginning and end of the transcript to get the effective length.
    Returns transcript tin score and its effective length.
    """

    # Get effective size
    # For short genes: if effective length <=0, then gene length is the effective length
    cvg = cvg[size:-size] if (size != 0 and len(cvg) - size * 2 > 0) else cvg

    eff_len = len(cvg)

    # Change to float
    cvg = list(map(float, cvg))

    # Remove zeros
    cvg = list(filter(None, cvg))

    if not cvg:
        tin = 0.0
        return tin, eff_len, 0

    # Calculate shannon's entropy
    ent = shannon_entropy(cvg)

    # Calculate uniformity
    uni = math.exp(ent)

    # Calculate tin on the effective length
    obs_tin = 100 * (uni) / eff_len

    # Calculate expected coverage
    exp_cov = sum(cvg) / eff_len

    # Calculate expected tin
    exp_tin = expected_tin(exp_cov)

    return obs_tin, exp_tin, eff_len


def expected_tin(cov):
    return (1 - math.exp(-cov)) * 100


def get_primary_aln(bam):
    """
    Extract primary alignments from bam.
    """

    pbam = get_filename(bam) + "_primary.bam"
    pbam = os.path.join(TMP, pbam)
    cmd1 = f'samtools view -h -b -F 4 -F 2304 {bam} >{pbam}'
    cmd2 = f'samtools index {pbam}'

    d1 = subprocess.run(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
                        universal_newlines=True)
    d2 = subprocess.run(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
                        universal_newlines=True)

    # Exit on error.
    exit_on_error(d1.stderr) if d1.stderr and d1.returncode != 0 else None
    exit_on_error(d2.stderr) if d2.stderr and d2.returncode != 0 else None
    return pbam


def store_background(fname):
    store = dict()

    stream = csv.reader(open(fname), delimiter="\t")
    for row in stream:
        tid, bg = row[:]
        store[tid] = bg
    return store


def subtract_background(elms, bg):
    def z(x):
        return 0 if x < 0 else x

    # Subtract background
    elms = [val - bg for val in elms]

    #  Make coverage to 0 if value < 0
    elms = list(map(z, elms))
    return elms


def get_effective_length(merged, size=50):
    """
    Get length as the sum of feature lengths.
    """

    flen = dict()

    def length(a):
        return (a[1] - a[0] + 1)

    for key, vals in merged.items():
        # size*2 since we subtract size bases from both ends while calculating tin.
        tot_len = sum(map(length, vals))
        eff_len = tot_len - size * 2
        eff_len = tot_len if eff_len <= 0 else eff_len
        flen[key] = eff_len
    # Header
    flen['samples'] = "eff_length"

    return flen


def exit_on_error(err):
    print(err)
    sys.exit()


def check_strand(strand):
    return True if strand in ['both', 'sense', 'antisense'] else False


def check_libtype(libtype):
    return True if libtype in ['single', 'paired'] else False


def check_bam(files):
    """
    Check if all files have the extension .bam
    """
    extns = [os.path.splitext(os.path.basename(f))[1] for f in files]
    return len(set(extns)) == 1


def get_extension(fname):
    """
    Returns first 3 characters in file extension.
    """
    ext = os.path.splitext(os.path.basename(fname))[1]
    # get 3 characters from extension.
    ext = ext[:4].lower()
    return ext


def check_ann(ann):
    """
    Check if annotation is a GTF/GFF/GFF3 file.
    """
    ext = get_extension(ann)
    if ext not in [".gtf", ".gff", ".gff3"]:
        return False
    return True


def check_inputs(bams, strand, lib_type, ann):
    if not (ann and bams):
        print("Alignment file(s) and annotation file in GTF/GFF3 format are required.")
        sys.exit()

    # Make sure the alignment file extension is bam.
    if not check_bam(bams):
        print("All alignment files must be have '.bam' extension.")
        sys.exit()

    if not check_ann(ann):
        print("GTF/GFF annotation file with gene feature row is required.")
        sys.exit()

    # Make sure entered strand is valid.
    if not check_strand(strand):
        print("Invalid strand. Available options are both, sense, antisense")
        sys.exit()

    # Make sure entered libtype is valid.
    if not check_libtype(lib_type):
        print("Invalid lib_type. Available options are single or paired")
        sys.exit()

    return


def get_random_string(length):
    # Choose from all lowercase letter.
    letters = string.ascii_lowercase
    random_str = ''.join(random.choice(letters) for i in range(length))
    return random_str


@plac.pos('bams', "comma separated bam files.")
@plac.opt('ann', help="annotation file in GTF or GFF3 format")
@plac.opt('feat', type=str, help="feature in the 3rd column of the anntation file on which TIN needs to be calculated")
@plac.opt('groupby', type=str, help="attribute by which features need to be combined , eg: gene_id")
@plac.opt('strand', type=str, help="strand on which tin should be calculated; options are both, sense or antisense")
@plac.opt('libtype', type=str, help="library type; options are paired or single")
@plac.flg('bg', help="when specified background noise will be subtracted")
@plac.opt('n', type=int, help="""number of bases to be subtracted from
                              each ends of the feature to calculate effective length""")
def run(bams, ann="", feat='exon', groupby='gene_id', strand='both', libtype='single', bg=False, n=50):
    bg_file, intron_len, intron_gtf = None, None, None

    # Check if inputs are valid.
    check_inputs(bams, strand, libtype, ann)

    # Create TMP directory
    os.makedirs(TMP)

    # Create run-log.
    log = os.path.join(TMP, "runlog.txt")

    # Make a string of bam files.
    bams = bams.split(',')
    bam = " ".join(bams)

    # Extract features from annotation file.
    genes, features = collect_features(ann=ann, feat_type=feat, groupby=groupby)

    # Merge overlapping features by groupby attribute.
    merged = merge_features(features)

    # Make a bed file of merged features.
    bed = make_feature_bed(genes, merged)

    # Get feature length
    feat_len = get_effective_length(merged, n)

    # Get intron coordinates for background noise calculation.
    if bg:
        introns = get_intron_coords(merged)
        intron_len = get_effective_length(introns, size=0)
        intron_gtf = make_intron_gtf(genes, introns)

    # Paired or single.
    paired = True if libtype == "paired" else False

    # Get counts.
    gene_counts = get_counts(bam=bam, ann=ann, strand=strand,
                             feat=feat, groupby=groupby, paired=paired, log=log)

    # Sample specific TIN calculations begin.
    obs_tins, exp_tins = dict(), dict()

    for idx, bam in enumerate(bams):

        # Calculate background noise if "--bg" is specified.
        if bg:
            bg_file = get_background_noise(bam=bam, strand=strand,
                                           intron_ann=intron_gtf, intron_len=intron_len, groupby=groupby, paired=paired)

        # Extract primary alignments
        pbam = get_primary_aln(bam)

        gene_tin = get_gene_tin(bam=pbam, bed=bed, strand=strand, bgfile=bg_file, size=n)

        # Collect obs_tins.
        for uid, vals in gene_tin.items():
            if uid == "samples":
                obs_tins.setdefault(uid, []).append(vals[0])
                exp_tins.setdefault(uid, []).append(vals[1])
            else:
                obs_tins.setdefault(uid, []).append(vals[2])
                exp_tins.setdefault(uid, []).append(vals[3])

    # Collect all results, gid,counts, exp_tin, obs_tin
    results = dict()

    for k, v in gene_counts.items():
        vals = list()
        vals.append(feat_len.get(k, ""))
        vals.extend(v)
        vals.extend(exp_tins.get(k, ""))
        vals.extend(obs_tins.get(k, ""))
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
    os.system(cmd)


if __name__ == "__main__":
    # plac.call(run)
    run()
