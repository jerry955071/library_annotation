"""
Design Sanger sequencing primer for gapped fastq using Primer3
"""
from io import TextIOWrapper
from Align.AlignIO import AlignIO
from Align.AlignRecord import AlignRecord
from typing import List
from pathlib import Path
import time
import argparse
import subprocess
import sys
import os


# argparse
parser = argparse.ArgumentParser()
# Primer3
parser.add_argument(
    "--p3",
    help="path to Primer3 executable",
    type=str
)
parser.add_argument(
    "--p3_settings_file",
    help="path to Primer3 setting file",
    type=str
)
# INPUTS
parser.add_argument(
    "--gapped_sam",
    metavar="GAPPED_SAM",
    help="filtered alignment results of GAPPED_FASTQ to reference CDS_FASTA"
)
# OUTPUTS
parser.add_argument(
    "--out_dir",
    help="file to write Primer3 output to",
    type=str
)
# PARAMS
parser.add_argument(
    "--palign_cutoff",
    help="minimal aligned percentages",
    type=float,
    default=0.49
)
parser.add_argument(
    "--len_noisy",
    help="",
    type=int,
    default=100
)
parser.add_argument(
    "--min_ovlp",
    help="minimal desired overlap between sequences",
    type=int,
    default=30
)
parser.add_argument(
    "--min_pickable",
    help="minimal sequences length used to pick primer",
    type=int,
    default=30
)


class primer3_core:
    def __init__(self, p3core:str, p3config:str) -> None:
        super().__init__()
        self.p3core = p3core
        self.p3config = p3config

    
    def run(self, records:List[dict]) -> List[dict]:
        tstamp = int(time.time()*1000)
        tmpin = f".{tstamp}_primer3_input.tmp"
        tmpout = f".{tstamp}_primer3_output.tmp"
        # input to p3
        with open(tmpin, "w") as handle:
            BoulderIO.write(handle, records)

        # run p3
        subprocess.run([
            self.p3core,
            f"--p3_settings_file={self.p3config}",
            f"--output={tmpout}",
            f"{tmpin}"
        ])

        # read p3 output
        with open(tmpout, "r") as handle:
            ret = BoulderIO.parse(handle)

        # remove tmp files
        os.remove(tmpin)
        os.remove(tmpout)

        return ret



class BoulderIO:    
    @staticmethod
    def parse(handle:TextIOWrapper) -> List[dict]:
        records = []
        record = {}
        for line in handle:
            line = line.strip()
            if line == "=":
                records.append(record)
                record = {}
            else:
                k, v = line.split("=")
                record[k] = v

        return records
    
    @staticmethod
    def write(handle:TextIOWrapper, records:List[dict]|dict):
        if isinstance(records, dict):
            for k, v in records.items():
                handle.write(f"{k}={v}\n")
            handle.write("=\n")
            return
        if isinstance(records, list):    
            for record in records:
                for k, v in record.items():
                    handle.write(f"{k}={v}\n")
                handle.write("=\n")
            return



# get/set arguments
args = parser.parse_args()
p3_core = primer3_core(args.p3, args.p3_settings_file)
gapped_sam = args.gapped_sam
out_dir = args.out_dir
palign_cutoff = args.palign_cutoff
len_noisy = args.len_noisy
min_ovlp = args.min_ovlp
min_pickable = args.min_pickable

# makedir
os.makedirs(out_dir, exist_ok=True)


# parsing sam file
align_obj = AlignIO(gapped_sam)
rlen_dict = {} # Gene ID -> Max transcript length
align_dict = {} # Sample name -> [AlignRecord1,...]
for l in align_obj.AlignGenerator(skip_header=False):
    if not isinstance(l, AlignRecord):
        if l[:3] == "@SQ":
            rname, rlen = l.rstrip().split("\t")[1:3]
            gene_id = rname.split("|")[1]
            try:
                rlen_dict[gene_id] = max([int(rlen[3:]), rlen_dict[gene_id]])
            except KeyError:
                rlen_dict[gene_id] = int(rlen[3:])
    else:
        sample = l.qname[:6]
        if sample in align_dict.keys():
            align_dict[sample].append(l)
        else:
            align_dict[sample] = [l]



# for align records in `align_dict`
for sample, aligns in align_dict.items():
    # Checking the goodness of align records
    # 1. Must contain just 2 records
    # 2. Both alignment results must be good enough, meaning that:
    #    (1) Must be mapped 
    #    (2) Nicely mapped, meaning that max(qp_align, rp_align) > palign_cutoff
    #        qp_align: align_col/qlen
    #        rp_align: align_col/rlen 
    
    # 1. Must contain just 2 records
    if len(aligns) != 2:
        # report and continue
        print(f"{sample}: {len(aligns)} record", file=sys.stderr)
        continue

    # 2. Both alignment results must be good enough
    failed_any = False
    for align in aligns:
        # (1) Must be mapped
        if align.rname == "*":
            failed_any = True
            print(f"{sample}: unmapped record {align.qname}", file=sys.stderr)
            continue

        # (2) Nicely mapped
        align_col = align.get_align_bln()
        qp_align = align_col / len(align.seq)
        rp_align = align_col / rlen_dict[align.rname.split("|")[1]]
        if not max(qp_align, rp_align) > palign_cutoff:
            failed_any = True
            print(f"{sample}: low mapping quality record {align.qname}", file=sys.stderr)
            continue
    # failed any
    if failed_any:
        continue

    # Both alignments are good
    fout = Path(f"{out_dir}/{sample}.primer3")
    open(fout, "w").close()
    front, back = aligns
    # n-subtituted sequence
    fn_sub = front.n_substitued()
    bn_sub = back.n_substitued()
    # first un-clipped base on query
    fqstart = front.get_softclip()[0]
    bqstart = back.get_softclip()[0]

    # 1. Get expected CDS length, identify gap on reference
    # A. if mapped to different reference CDS
    if front.rname != back.rname:
        fnd = len(fn_sub) - min_ovlp - len_noisy
        fst = fnd - 200
        bst = 0 + min_ovlp + len_noisy
        bnd = bst + 200
        expected_gap = \
            rlen_dict[front.rname.split("|")[1]] - \
                len(front.seq) + sum(front.get_softclip()) - \
                    len(back.seq) + sum(back.get_softclip())
        # Primer3
        lp3record = p3_core.run({
            "SEQUENCE_ID": f"{front.qname}[{fqstart + fst}:{fqstart + fnd}]",
            "SEQUENCE_TEMPLATE": fn_sub[fst:fnd],
            "PRIMER_PICK_LEFT_PRIMER": 1
        })[0]
        rp3record = p3_core.run({
            "SEQUENCE_ID": \
                f"{back.qname}[{bqstart + bst}:{bqstart + bnd}]",
            "SEQUENCE_TEMPLATE": bn_sub[bst:bnd],
            "PRIMER_PICK_RIGHT_PRIMER": 1
        })[0]
        # A1
        if expected_gap > 1050 - 2 * (min_ovlp + len_noisy):
            # write both primer to file
            with open(fout, "a") as handle:
                BoulderIO.write(handle, lp3record)
                BoulderIO.write(handle, rp3record)
            print(f"{sample}: diff ref, 2D strategy,  expected gap = {expected_gap} bp", file=sys.stdout)
            continue
        else:
        # A2
            # pick the better primer, if exist go next loop
            D = {}
            try:
                D[float(lp3record["PRIMER_LEFT_0_PENALTY"])] = lp3record
            except KeyError:
                D[9999] = lp3record
            try:
                D[float(rp3record["PRIMER_RIGHT_0_PENALTY"])] = rp3record
            except KeyError:
                D[9999] = rp3record
        
            with open(fout, "a") as handle:
                BoulderIO.write(handle, D[min(D.keys())])
            
            print(f"{sample}: diff ref, 1D strategy, expected gap = {expected_gap} bp", file=sys.stdout)
            continue
            
    # if mapped to the same reference CDS
    expt_cds_len = rlen_dict[front.rname.split("|")[1]]
    # gap_start = map end front
    gap_start = front.pos + front.cigar.count("=XD")
    # gap_end = map start back
    gap_end = back.pos
    gap_mid = int((gap_start + gap_end) / 2)
    gap_on_ref = gap_end - gap_start
    # NOTE: SPECIAL CASES
    if gap_on_ref <= 0:
        print(f"{sample}: {abs(gap_on_ref)} bp overlap on the same ref", file=sys.stderr)
        os.remove(fout)
        continue
    
    # 2. Design primer using Primer3
    fnd = gap_start - min_ovlp - len_noisy
    bst = gap_end + min_ovlp + len_noisy
    # Try 1-way strategy 
    if fnd - bst + 1050 >= min_pickable:
        fst = max(bst - 1050, 0)
        bnd = fnd + 1050 # might exceed read length 
        # Primer3
        lp3record = p3_core.run({
            "SEQUENCE_ID": f"{front.qname}[{fqstart + fst}:{fqstart + fnd}]",
            "SEQUENCE_TEMPLATE": fn_sub[fst - front.pos + 1:fnd - front.pos + 1],
            "PRIMER_PICK_LEFT_PRIMER": 1
        })[0]
        rp3record = p3_core.run({
            "SEQUENCE_ID": \
                f"{back.qname}[{bqstart + bst - gap_end}:{bqstart + bnd - gap_end}]",
            "SEQUENCE_TEMPLATE": bn_sub[bst - back.pos + 1:bnd - back.pos + 1],
            "PRIMER_PICK_RIGHT_PRIMER": 1
        })[0]
        
        # pick the better primer, if exist go next loop
        D = {}
        try:
            D[float(lp3record["PRIMER_LEFT_0_PENALTY"])] = lp3record
        except KeyError:
            D[9999] = lp3record
        try:
            D[float(rp3record["PRIMER_RIGHT_0_PENALTY"])] = rp3record
        except KeyError:
            D[9999] = rp3record
        
        if not min([D.keys()]) == 9999:
            # write to file
            with open(fout, "a") as handle:
                BoulderIO.write(handle, D[min(D.keys())])
            print(f"{sample}: same ref, 1D strategy, gap on ref = {gap_on_ref} bp", file=sys.stdout)
            continue # next loop
    
    # Try 2-way strategy
    if fnd - gap_mid - min_ovlp + 1050 >= min_pickable:
        lmid = rmid = True
        # Try 2-way mid-point strategy
        fst = gap_mid + min_ovlp - 1050
        bnd = gap_mid - min_ovlp + 1050
        # Primer3
        lp3record = p3_core.run({
            "SEQUENCE_ID": f"{front.qname}[{fqstart + fst}:{fqstart + fnd}]",
            "SEQUENCE_TEMPLATE": fn_sub[fst - front.pos + 1:fnd - front.pos + 1],
            "PRIMER_PICK_LEFT_PRIMER": 1
        })[0]
        rp3record = p3_core.run({
            "SEQUENCE_ID": \
                f"{back.qname}[{bqstart + bst - gap_end}:{bqstart + bnd - gap_end}]",
            "SEQUENCE_TEMPLATE": bn_sub[bst - back.pos + 1:bnd - back.pos + 1],
            "PRIMER_PICK_RIGHT_PRIMER": 1
        })[0]

        # If no primer returned, extend region to 200 bp
        if lp3record["PRIMER_LEFT_NUM_RETURNED"] == "0":
            if fst > fnd - 200:
                fst = fnd - 200
                lp3record = p3_core.run({
                    "SEQUENCE_ID": f"{front.qname}[{fqstart + fst}:{fqstart + fnd}]",
                    "SEQUENCE_TEMPLATE": fn_sub[fst - front.pos + 1:fnd - front.pos + 1],
                    "PRIMER_PICK_LEFT_PRIMER": 1
                })[0]
                lmid = False

        if rp3record["PRIMER_RIGHT_NUM_RETURNED"] == "0":
            if bnd < bst + 200:
                bnd = bst + 200
                rp3record = p3_core.run({
                    "SEQUENCE_ID": \
                        f"{back.qname}[{bqstart + bst - gap_end}:{bqstart + bnd - gap_end}]",
                    "SEQUENCE_TEMPLATE": bn_sub[bst - back.pos + 1:bnd - back.pos + 1],
                    "PRIMER_PICK_RIGHT_PRIMER": 1
                })[0]
                rmid = False

        # write to file
        with open(fout, "a") as handle:
            BoulderIO.write(handle, lp3record)
            BoulderIO.write(handle, rp3record)
        
        lstrat = "mid-point" if lmid else 200
        rstrat = "mid-point" if rmid else 200
        print(f"{sample}: same ref, 2D ({lstrat}, {rstrat}) strategy, gap on ref = {gap_on_ref} bp", file=sys.stdout)
        continue

    print(f"{sample}: no primer found", file=sys.stdout)
