"""
Design Sanger sequencing primer for gapped fastq using Primer3
"""
from Align.AlignIO import AlignIO
from Align.AlignRecord import AlignRecord
from Bio import SeqIO
import argparse
import subprocess
import sys
import os


# argparse
parser = argparse.ArgumentParser()
# Primer3
parser.add_argument(
    "--p3",
    help="path to Primer3 executable (primer3_core)",
    type=str
)
parser.add_argument(
    "--p3_settings_file",
    help="path to Primer3 setting file",
    type=str
)
# INPUTS
parser.add_argument(
    "--gapped_fq",
    help="directory containing gapped fastq",
    type=str
)
parser.add_argument(
    "--filtered_sam",
    metavar="FILTERED_SAM",
    help="filtered alignment results of RAW_FASTQ to reference CDS_FASTA"
)
# OUTPUTS
parser.add_argument(
    "--out_dir",
    help="file to write Primer3 output to",
    type=str
)
# PARAMS
parser.add_argument(
    "--offset",
    help="",
    type=int,
    default=150
)
parser.add_argument(
    "--max_interval",
    help="",
    type=int,
    default=200
)

_intervals = [60, 120, 160, 200]

# get arguments
args = parser.parse_args()
p3 = args.p3
p3_settings_file = args.p3_settings_file
gapped_fq = args.gapped_fq
filtered_sam = args.filtered_sam
out_dir = args.out_dir; os.makedirs(out_dir, exist_ok=True)
offset = args.offset
max_interval = args.max_interval


# parsing sam file
align_obj = AlignIO(filtered_sam)
rlen_dict = {} # Reference ID -> Reference length
align_dict = {} # Query name -> Align Record
for l in align_obj.AlignGenerator(skip_header=False):
    if not isinstance(l, AlignRecord):
        # rlen_dict
        if l[:3] == "@SQ":
            rname, rlen = l.rstrip().split("\t")[1:3]
            rlen_dict[rname[3:]] = int(rlen[3:])
    else:
        # align_dict
        align_dict[l.qname] = l


# for records in '../outputs/gapped_fastq'
for fname in sorted(os.listdir(gapped_fq)):
    # get base name
    bname = fname.strip(".fq")
    
    # read records
    with open(f"{gapped_fq}/{fname}", "r") as handle_in:
        records = [i for i in SeqIO.parse(handle_in, format="fastq")]
        
    ## NOTE: 'continue' if not 2 records
    if len(records) != 2:
        print(f"{fname}: {len(records)} records", file=sys.stderr)
        continue

    # get expected transcripts length
    ## NOTE: 'continue' if unmapped
    expt_cds = {align_dict[r.id].rname: r.id for r in records}
    if "*" in expt_cds.keys():
        tmp = expt_cds["*"]
        print(f"{fname}: unmapped record {tmp}", file=sys.stderr)
        continue
    expt_cds_len = max([rlen_dict[i] for i in expt_cds])

    # get expected gap on reference
    ## NOTE: 'continue' if overlapped on reference 
    aligns = [align_dict[record.id] for record in records]
    gap_start = aligns[0].pos + aligns[0].cigar.count("=XD")
    gap_end = aligns[1].pos
    gap_on_ref = gap_end - gap_start
    if gap_on_ref <= 0:
        print(f"{fname}: overlapped on reference", file=sys.stdout)
        continue

    # design primer using Primer3
    fwd_align, rev_align = aligns
    
    ## design left primer on forward sequence
    p3in = ".p3in.tmp"
    p3out = ".p3out.tmp"    
    n_subed = fwd_align.n_substitued()
    end = -offset
    bclip = 0 if fwd_align.cigar.cigar_str[-1] != "S" else \
        fwd_align.cigar.cigar_int[-1]
    x = len(fwd_align.seq) - bclip
    end_loop = False
    for interval in _intervals:
        start = end - interval
        with open(p3in, "w") as handle:
            handle.write(
                f"SEQUENCE_ID={fwd_align.qname}[{x + start}:{x + end}]\n"
                f"SEQUENCE_TEMPLATE={n_subed[start:end]}\n"
                "PRIMER_TASK=generic\n"
                "PRIMER_PICK_LEFT_PRIMER=1\n"
                "PRIMER_PICK_INTERNAL_OLIGO=0\n"
                "PRIMER_PICK_RIGHT_PRIMER=0\n"
                "PRIMER_PRODUCT_OPT_SIZE=30\n"
                "PRIMER_EXPLAIN_FLAG=1\n"
                "=\n"
            )
        proc = subprocess.run([
            p3,
            f"--p3_settings_file={p3_settings_file}",
            f"--output={p3out}",
            p3in
            ]
        )
        with open(p3out, "r") as handle:
            p3record = {}
            for line in handle:
                line = line.strip()
                if line != "=":
                    k, v = line.split("=")
                    p3record[k] = v
        
        # NOTE: 'break' if find left primer
        if p3record["PRIMER_LEFT_NUM_RETURNED"] != "0":
            break
    
    subprocess.run([
        "cp",
        p3out,
        f"{out_dir}/{bname}.primer3"
        ]
    )
    
    ## design right primer on reverse sequence
    if gap_on_ref > (1050 - offset * 2):
        p3in = ".p3in.tmp"
        p3out = ".p3out.tmp"    
        n_subed = rev_align.n_substitued()
        start = offset
        fclip = 0 if rev_align.cigar.cigar_str[0] != "S" else \
            rev_align.cigar.cigar_int[0]
        x = fclip
        end_loop = False
        for interval in _intervals:
            end = start + interval
            # write p3 input
            with open(p3in, "w") as handle:
                handle.write(
                    f"SEQUENCE_ID={rev_align.qname}[{x + start}:{x + end}]\n"
                    f"SEQUENCE_TEMPLATE={n_subed[start:end]}\n"
                    "PRIMER_TASK=generic\n"
                    "PRIMER_PICK_LEFT_PRIMER=0\n"
                    "PRIMER_PICK_INTERNAL_OLIGO=0\n"
                    "PRIMER_PICK_RIGHT_PRIMER=1\n"
                    "PRIMER_PRODUCT_OPT_SIZE=30\n"
                    "PRIMER_EXPLAIN_FLAG=1\n"
                    "=\n"
                )
            # execute p3
            proc = subprocess.run([
                p3,
                f"--p3_settings_file={p3_settings_file}",
                f"--output={p3out}",
                p3in
                ]
            )
            # read p3 output
            with open(p3out, "r") as handle:
                p3record = {}
                for line in handle:
                    line = line.strip()
                    if line != "=":
                        k, v = line.split("=")
                        p3record[k] = v
                
            # NOTE: 'break' if find left primer
            if p3record["PRIMER_RIGHT_NUM_RETURNED"] != "0":
                break
        
        with open(f"{out_dir}/{bname}.primer3", "a") as hout, \
            open(p3out) as hin:
            for l in hin:
                hout.write(l)
        

os.remove(p3in)
os.remove(p3out)