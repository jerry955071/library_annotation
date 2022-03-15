from Bio import SeqIO
import argparse


# 
parser = argparse.ArgumentParser(description='Reminder:')
parser.add_argument(
    "-i",
    metavar="",
    nargs="+",
    help="input fastq files (separated by space)",
    type=str
    )
parser.add_argument(
    "--out_fa",
    metavar="",
    help="combined fasta records file",
    type=str
    )


#
args = parser.parse_args()
flist = args.i
out_fa = args.out_fa


#
open(out_fa, "w").close()
with open(out_fa, "a") as handle_out:
    for fname in flist:
        basename = fname.split("/")[-1].split(".")[0]
        with open(fname, "r") as handle_in:
            for record in SeqIO.parse(handle_in, "fastq"):
                record.description = basename
                SeqIO.write(record, handle_out, "fasta")

