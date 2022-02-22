"""
Extract CDS from Gencode reference (protein coding transcripts)
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import re
from typing import Tuple


# create ArgumentParser
parser = argparse.ArgumentParser(description='Reminder:')
parser.add_argument("-i", help="protein coding transcripts", type=str)
parser.add_argument("-o", help="output coding sequences", type=str)


# get parsed command line arguments
args = parser.parse_args()
fin = args.i
fout = args.o
# fin = "../genomic_data/gencode.v38.pc_transcripts.fa"
# fout = "../genomic_data/gencode.v38.cds.fa"


# create output file
open(fout, "w").close()


# get range of CDS from string
def where_cds(x:str) -> Tuple[int]:
    p_cds = re.compile("CDS:[0-9]+-[0-9]+")
    p_int = re.compile("[0-9]+")
    match_cds = re.findall(p_cds, x)[0]
    match_int = tuple(int(i) for i in re.findall(p_int, match_cds))
    return match_int


# extract cds
with open(fin, "r") as handle_in, open(fout, "a") as handle_out:
    for record in SeqIO.parse(handle_in, "fasta"):
        cds = where_cds(record.id)
        record.seq = record.seq[cds[0] - 1:cds[1]]
        foo = SeqIO.write(record, handle_out, "fasta") \
            # foo catches output `1`
