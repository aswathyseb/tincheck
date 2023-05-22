import plac
import argparse, subprocess, string


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


def write_to_file(rows, fh):
    for row in rows:
        fh.write(row + "\n")


@plac.pos('bams', "bam files")
@plac.opt('outfile', "Output file name.")
def run(bams, outfile='ann.gtf'):
    bam = bams.split(',')[0]
    refs = parse_bam_header(bam)
    start, strand = 1, "+"
    fh = open(outfile, "w")
    for r in refs:
        chrom, end = r[0], r[1]
        gene_ln = make_gtf_feature_rows(chrom, start, end, strand, chrom, idx=0, feat="gene")
        trans_ln = make_gtf_feature_rows(chrom, start, end, strand, chrom, idx=0, feat="transcript")
        exon_ln = make_gtf_feature_rows(chrom, start, end, strand, chrom, idx=0, feat="exon")
        write_to_file([gene_ln, trans_ln, exon_ln], fh)


if __name__ == "__main__":
    # plac.call(run)
    run()
