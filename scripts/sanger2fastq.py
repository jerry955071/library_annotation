"""
1. Convert AB1 to FASTQ.
2. Replace file name typo ('-' -> '_')
3. Reverse complement reverse sequences (recognized by *PATTERN*)
"""

from Bio import SeqIO
import argparse
import os
import sys
import re

# create ArgumentParser
parser = argparse.ArgumentParser(description='Reminder:')
parser.add_argument("--in_dir", help="input .ab1 directory", type=str)
parser.add_argument("--out_dir", help="output .fq directory", type=str)
parser.add_argument(
    "--pattern_rev",
    help="the pattern used to recognize reverse sequencese",
    type=str,
    default=".*(ADR).*"
    )

# get parsed command line arguments
args = parser.parse_args()
in_dir = args.in_dir
out_dir = args.out_dir
pattern_rev = args.pattern_rev

# create reverse file name parser
prog = re.compile(pattern_rev)

# list all sequencing results files
flist = os.listdir(in_dir)
flist.sort()

# create output directory (exist_ok)
os.makedirs(out_dir, exist_ok=True)

# 
n_abi = 0
n_fq = 0
for fname in flist:
    # skip not '.ab1' files
    if not ".ab1" in fname:
        continue
    
    n_abi += 1
    # get file base name
    bname = fname.split(".")[0]

    # 'abi' to 'fastq'
    with open(f"{out_dir}/{bname}.fq", "w") as handle_out:
        # get abi record
        record = SeqIO.read(f"{in_dir}/{fname}", "abi")
        
        # replace typo
        if "-" in record.name:
            record.id = record.name.replace("-", "_")
        else:
            record.id = record.name

        # reverse reverse sequences
        if prog.match(record.id):
            rev = record.reverse_complement()
            rev.id = record.id
            rev.name = record.name
            rev.description = ""
            record = rev
        
        # output to 'fastq'
        foo = SeqIO.write(record, handle_out, "fastq")
        n_fq += 1

print(f"""\
INPUT DIRECTORY:\t'{in_dir}'
OUTPUT DIRECTORY:\t'{out_dir}'
INPUT AB1 FILES:\t{n_abi}
OUTPUT FASTQ FILES:\t{n_fq}""", file=sys.stdout
    )
        
