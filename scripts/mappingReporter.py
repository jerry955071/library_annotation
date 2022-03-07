from .Align.AlignIO import AlignIO
import argparse
import re
import pandas as pd


# create ArgumentParser
parser = argparse.ArgumentParser(description='Reminder:')
parser.add_argument(
    "-i",
    metavar="",
    help="alignment result .sam file",
    type=str
    )
parser.add_argument(
    "--out_summary", 
    metavar="",
    help="output summary file (.tsv)", 
    type=str
    )
parser.add_argument(
    "--out_hist",
    metavar="",
    help="output histogram in text format (.txt)", 
    type=str
    )
parser.add_argument(
    "--out_unmapped",
    metavar="",
    help="output unmapped read names (one read per line)", 
    type=str
    )


# get parsed command line arguments
args = parser.parse_args()
fin = args.i
out_summary = args.out_summary
out_hist = args.out_hist
out_unmapped = args.out_unmapped


# get range of CDS from Gencode format GeneID
def get_cds_len(x:str):
    p_cds = re.compile("CDS:[0-9]+-[0-9]+")
    p_int = re.compile("[0-9]+")
    match_cds = re.findall(p_cds, x)[0]
    match_int = tuple(int(i) for i in re.findall(p_int, match_cds))
    return match_int[1] - match_int[0] + 1


# initiate summary_dict
summary_dict = {
    "qname": [],
    "cds_start": [],
    "cds_end": [],
    "map_start": [],
    "map_end": [],
    "pid": [],
    "pid_adj": [],
    "n_error": [],
    "n_insert": [],
    "n_mismatch": [],
    "n_deletion": [],
    "rname": []
}


# initiate unmapped.fq
open(out_unmapped, "w").close


with open(fin, "r") as handle:
    for record in AlignIO(handle).AlignGenerator():
        # if unmapped
        if record.rname == "*":
            # output unmapped records to file
            with open(out_unmapped, "a") as fout:
                fout.write(f"{record.qname}\n")
                fout.write(f"{record.seq}\n")
                fout.write("+\n")
                fout.write(f"{record.qual}\n")

        else:
            # get 
            cds_len = get_cds_len(record.rname)
            cds_reg = (1, cds_len)
            map_col = record.calculate("map_col")
            n_err = record.calculate("n_error")
            # NOTE: formula of percent identity (PID) and adj. PID 
            pid = 100 * record.cigar.count("=") / record.cigar.count("=XID")
            pid_adj = pid * record.cigar.count("=XD") / cds_len
            
            # catch duplicated record
            if record.qname in summary_dict["qname"]:
                raise Exception(f"Found duplicated record: {record.qname}")

            # expand `summary_dict`
            summary_dict["qname"].append(record.qname)
            summary_dict["rname"].append(record.rname)
            summary_dict["cds_start"].append(cds_reg[0])
            summary_dict["cds_end"].append(cds_reg[1])
            summary_dict["map_start"].append(map_col[0])
            summary_dict["map_end"].append(map_col[1])
            summary_dict["pid"].append("%.1f" % pid)
            summary_dict["pid_adj"].append("%.1f" % pid_adj)
            summary_dict["n_error"].append(n_err)
            summary_dict["n_insert"].append(record.cigar.count("I"))
            summary_dict["n_mismatch"].append(record.cigar.count("X"))
            summary_dict["n_deletion"].append(record.cigar.count("D"))


# output `summary_dict` to file
df = pd.DataFrame(summary_dict)
df.to_csv(out_summary, sep="\t", index=False)


# create histogram
def hist(
    x:list, 
    interval:tuple, 
    binwidth:float,
    title:str="",
    outfile:str=None
    ) -> None:
    
    start, end = interval
    # initiate `hist_dict`
    hist_dict = {}
    k = start
    while k <= end:
        left = k
        right = k + binwidth
        key = f"{k:.1f}"
        hist_dict[key] = 0
        k += binwidth

    for i in x:
        bins = f"{float(i) // binwidth * binwidth:.1f}"
        hist_dict[bins] += 1
    
    # output
    if outfile == None:
        print(f"{title}")
        print(f"Histogram:")
        k = start
        while k <= end:
            left = k
            right = k + binwidth
            right = right if right <= end else left
            key = f"{k:.1f}"
            print(f"\t{left:.1f}-{right:.1f}: {hist_dict[key]}")
            k += binwidth
        print("\nParams:")
        print(f"\t-interval: ({start}, {end})")
        print(f"\t-binwidth: {binwidth}\n")
    else:
        with open(outfile, "a") as fout:
            fout.write(f"{title}\n")
            fout.write(f"Histogram:\n")
            k = start
            while k <= end:
                left = k
                right = k + binwidth
                right = right if right <= end else left
                key = f"{k:.1f}"
                fout.write(f"\t{left:.1f}-{right:.1f}: {hist_dict[key]}\n")
                k += binwidth
            fout.write("\nParams:\n")
            fout.write(f"\t-interval: ({start}, {end})\n")
            fout.write(f"\t-binwidth: {binwidth}\n\n")
    return


# initiate `out_hist`
open(out_hist, "w").close

# append histogram to file
hist(
    df["pid"],
    interval=(0, 100),
    binwidth=5.0,
    outfile=out_hist,
    title="Histogram of percent identity"
    )

hist(
    df["pid_adj"],
    interval=(0, 100),
    binwidth=5.0,
    outfile=out_hist,
    title="Histogram of adjusted percent identity"
    )

