from Bio import SeqIO
import argparse
import os

# create ArgumentParser
parser = argparse.ArgumentParser(description='Reminder:')
parser.add_argument("-i", 
    help="folder containing sanger sequencing results", type=str)
parser.add_argument("--out1", help="output read1 file name", type=str)
parser.add_argument("--out2", help="output read2 file name", type=str)
parser.add_argument("--suff1", help="input read1 file suffix", type=str)
parser.add_argument("--suff2", help="input read2 file suffix", type=str)


# get parsed command line arguments
args = parser.parse_args()
in_dir = args.i
out1 = args.out1
out2 = args.out2
suff1 = args.suff1
suff2 = args.suff2
# in_dir = "raw_sanger"
# out1 = "raw_fastq/read1.fq"
# out2 = "raw_fastq/read2.fq"
# suff1 = "ADF.ab1"
# suff2 = "ADR.ab1"


# create files
handle1 = open(out1, "w").close()
handle2 = open(out2, "w").close()


# check sequencing results files
flist = os.listdir(in_dir)
flist.sort()


# write .ab1 to .fq
with open(out1, "a") as handle1, open(out2, "a") as handle2:
    flist = os.listdir(in_dir)
    flist.sort()
    for fname in flist:
        if suff1 in fname:
            record = SeqIO.read(f"{in_dir}/{fname}", "abi")
            # some record.id were mis-named use record.name instead
            record.id = record.name
            # remove suffix
            record.id = record.id.split("_ADF")[0]
            # replace "-" with "_" in mis-named file names
            if "-" in record.id:
                record.id = record.id.replace("-", "_")
            # output to file
            foo = SeqIO.write(record, handle1, "fastq") \
                # `foo` were used to catch the output `1` from `SeqIO.write()`

        if suff2 in fname:
            record = SeqIO.read(f"{in_dir}/{fname}", "abi")
            # some record.id were mis-named use record.name instead
            record.id = record.name
            # remove suffix
            record.id = record.id.split("_ADR")[0]
            # replace "-" with "_" in mis-named file names
            if "-" in record.id:
                record.id = record.id.replace("-", "_")
            # output to file
            foo = SeqIO.write(record, handle2, "fastq") \
                # `foo` were used to catch the output `1` from `SeqIO.write()`
