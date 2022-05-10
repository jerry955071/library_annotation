from Bio import SeqIO
import mergeSeqRecord
import sys
import os
import re
import argparse
from typing import List


# create ArgumentParser
parser = argparse.ArgumentParser(description='Parameters for this script') 
parser.add_argument(
    "-w",
    help="window size of the minimizer algorithm for detecting overlap",
    default=13,
    type=int
)
parser.add_argument(
    "-k",
    help="kmer size of the minimizer algorithm for detecting overlap",
    default=8,
    type=int
)
parser.add_argument(
    "--min_ovlp",
    help="minimal required overlap length",
    default=30,
    type=int
)
parser.add_argument(
    "--in_dir",
    help="input fastq directory",
    type=str,
    default=None
)
parser.add_argument(
    "--merged_dir",
    help="output merged fastq directory",
    type=str,
    default=None
)
parser.add_argument(
    "--gapped_dir",
    help="output gapped fastq directory",
    type=str,
    default=None
)
parser.add_argument(
    "--report",
    help="output report file name",
    type=str,
    default=None
)


# file name pattern
pattern = \
    "(?P<sample>[0-9]{2}_[A-Z][0-9]{2})_(?P<infix>.*)[.](?P<suff>f(ast)?q)"
prog = re.compile(pattern)

# file name sorter
def sort_name_list(fnames:List[str]) -> List[str]:
    # pattern = \
    # "(?P<sample>[0-9]{2}_[A-Z][0-9]{2})_(?P<infix>.*)[.](?P<suff>f(ast)?q)"
    # prog = re.compile(pattern)
    # fnames.sort(key=lambda x:prog.match(x).group("infix"))
    return fnames.sort()


# get arguments from argparse 
args = parser.parse_args()
w = args.w
k = args.k
min_ovlp = args.min_ovlp
in_dir=args.in_dir
merged_dir=args.merged_dir
gapped_dir=args.gapped_dir
out_report = args.report if args.report else sys.stderr

# create output directory
os.makedirs(merged_dir, exist_ok=True)
os.makedirs(gapped_dir, exist_ok=True)

# get list of input fq
flist = os.listdir(in_dir)
flist.sort()

# create file map
file_map: dict[str, list]
file_map = {}
for f in flist:
    result = prog.match(f)
    if not result:
        continue
    
    sample = result.group("sample")
    if sample in file_map.keys():
        file_map[sample].append(f)
    else:
        file_map[sample] = [f] 


# for each key-value pairs in `file_map`
for sample in file_map.keys():
    # sort file names
    file_map[sample] = sort_name_list(file_map[sample])

    # get records
    records = []
    for f in file_map[sample]:
        record = SeqIO.read(f"{in_dir}/{f}", format="fastq")
        records.append(record)

    # only 1 record
    if len(records) == 1:
        with open(f"{gapped_dir}/{sample}.fq", "w") as handle:
            SeqIO.write(records[0], handle, "fastq")
    
    # more than 1 records
    else:
        # merge records
        i = 0
        while i < len(records) - 1:
            # get fastq
            fq1 = records[i]
            fq2 = records[i + 1]

            # try merge 2 fastq
            merged = mergeSeqRecord.tryMerge(fq1, fq2, w, k, min_ovlp)
            
            # if merged
            if merged:
                records.pop(i)
                records.pop(i)
                records.insert(i, merged)
            else:
                i += 1

        # if all merged, output to 'merged_dir'
        if len(records) == 1:
            with open(f"{merged_dir}/{sample}.fq", "w") as handle:
                SeqIO.write(records[0], handle, "fastq")
        
        # else output to 'gapped_dir'
        else:
            open(f"{gapped_dir}/{sample}.fq", "w").close()
            with open(f"{gapped_dir}/{sample}.fq", "a") as handle:    
                for record in records:
                    foo = SeqIO.write(record, handle, "fastq")

