from Bio import SeqIO
import mergeSeqRecord
import argparse
import pandas as pd


# create ArgumentParser
parser = argparse.ArgumentParser(description='Parameters for this script') 
parser.add_argument(
    "-w",
    help="window size for detecting overlap",
    default=13,
    type=int
)
parser.add_argument(
    "-k",
    help="kmer size for detecting overlap",
    default=8,
    type=int
)
parser.add_argument(
    "--min_ovlp",
    help="minimal required overlap length",
    default=30,
    type=int
)
parser.add_argument("--out_merged", help="merged fastq filename", type=str)
parser.add_argument("--out_single", help="single-end fastq filename", type=str)
parser.add_argument("--out_paired1", help="paired-end read1 filename", type=str)
parser.add_argument("--out_paired2", help="paired-end read2 filename", type=str)
parser.add_argument("--in1", help="input read1 file name", type=str)
parser.add_argument("--in2", help="input read2 file name", type=str)
parser.add_argument("--summary", help="output summary file name", type=str)
parser.add_argument("--report", help="output report file name", type=str)


# get arguments from argparse 
args = parser.parse_args()
w =args.w
k = args.k
min_ovlp = args.min_ovlp
out_merged = args.out_merged
out_single = args.out_single
out_paired1 = args.out_paired1
out_paired2 = args.out_paired2
in1 = args.in1
in2 = args.in2
out_report = args.report
out_summary = args.summary

# # for test
# out_merged = "merged.tmp"
# out_single = "single.tmp"
# out_paired1 = "paired1.tmp"
# out_paired2 = "paired2.tmp"
# in1 = "../output_fastp/forward_passed.fq"
# in2 = "../output_fastp/reverse_passed.fq"
# out_report = "report.tmp"
# out_summary = "summary.tmp"


# access fastq files using `SeqIO.index()`
record_dict1 = SeqIO.index(in1, "fastq")
record_dict2 = SeqIO.index(in2, "fastq")


# create output files (clear file if already exists)
open(out_merged, "w").close()
open(out_single, "w").close()
open(out_paired1, "w").close()
open(out_paired2, "w").close()
open(out_report, "w").close()
open(out_summary, "w").close()


# initiate summary_dict
summary_dict = {
    "rname":[],
    "forward":[],
    "reverse":[],
    "overlapped":[]
}


# get all sample names
all_samples = set(record_dict1._offsets.keys()) | \
    set(record_dict2._offsets.keys())
all_samples = list(all_samples)
all_samples.sort()


# category counter
report_dict = {
    "se": int(),
    "pe": int(),
    "pe_ovlp": int(),
    "pe_pair": int()
}


# for each sample do:
for name in all_samples:
    if name == "01_B05":
        pass
    # check if forward/reverse sequencing results exists
    fwd_exist = name in set(record_dict1._offsets.keys())
    rev_exist = name in set(record_dict2._offsets.keys())

    # expand report_dict
    summary_dict["rname"].append(name)
    summary_dict["forward"].append(fwd_exist)
    summary_dict["reverse"].append(rev_exist)

    # if both forward/reverse sequencing results exist
    if fwd_exist and rev_exist:
        # record read type to `report_dict`
        report_dict["pe"] += 1

        # get forward/reverse sequences
        fwd = record_dict1[name]
        rev = record_dict2[name].reverse_complement()
        rev.id = record_dict2[name].id
        rev.description = ""


        # try merging two reads
        merged = mergeSeqRecord.tryMerge(fwd, rev, w, k, min_ovlp)
        
        if merged is None:
            # record read type to `report_dict`
            report_dict["pe_pair"] += 1

            # record to `summary_dict`
            summary_dict["overlapped"].append(False)

            # output to file
            with open(out_paired1, "a") as handle:
                foo = SeqIO.write(fwd, handle, "fastq")
            
            with open(out_paired2, "a") as handle:
                foo = SeqIO.write(rev, handle, "fastq")
            
            continue
 
        else:
            # record read type to `report_dict`
            report_dict["pe_ovlp"] += 1

            # record to `summary_dict`
            summary_dict["overlapped"].append(True)

            # output to file
            with open(out_merged, "a") as handle:
                foo = SeqIO.write(merged, handle, "fastq")
            
            continue

    # if only forward or revese sequencing result exist
    if fwd_exist:
        fwd = record_dict1[name]

        # record read type to `report_dict`
        report_dict["se"] += 1

        # record to `summary_dict`
        summary_dict["overlapped"].append(False)

        # output to file
        with open(out_single, "a") as handle:
                foo = SeqIO.write(fwd, handle, "fastq")
        
        continue
    
    if rev_exist:
        rev = record_dict2[name].reverse_complement()
        rev.id = record_dict2[name].id
        rev.description = ""

        # record read type to `report_dict`
        report_dict["se"] += 1

        # record to `summary_dict`
        summary_dict["overlapped"].append(False)

        # output to file
        with open(out_single, "a") as handle:
                foo = SeqIO.write(rev, handle, "fastq")

        continue

# output to summary.tsv
pd_summary = pd.DataFrame(summary_dict)
pd_summary.to_csv(out_summary, sep="\t", index=False)


# output report.txt
with open(out_report, "a") as handle:
    handle.write(
        f"Number of single-end (SE) reads: {report_dict['se']}\n" + 
        f"Number of paired-end (PE) reads: {report_dict['pe']}\n" + 
        f"         ─ overlapping PE reads: {report_dict['pe_ovlp']}\n" + 
        f"     ─ not-overlapping PE reads: {report_dict['pe_pair']}\n"
        )


# need_re_seq = []
# for index, row in pd_report.iterrows():
#     if row["forward"] is True and row["reverse"] is True:
#         if row["overlapped"] is False:
#             need_re_seq.append(row["rname"])
