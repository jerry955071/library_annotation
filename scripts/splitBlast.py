import argparse
from io import TextIOWrapper


# create `ArgumentParser`
parser = argparse.ArgumentParser(description='Parameters for this script')
parser.add_argument("-i", help="input blast output (fmt = 0)", type=str)
parser.add_argument("--out_dir", help="output directory", type=str)


# get arguments from argparse 
args = parser.parse_args()
i =args.i
out_dir = args.out_dir


out_handle = 0
with open(i) as fin:
    for line in fin:
        if "Query= " in line:
            fname = line.strip().split(" ")[1] + ".txt"
            out_handle = open(f"{out_dir}/{fname}", "w")
            out_handle = open(f"{out_dir}/{fname}", "a")
            out_handle.write(line)
        elif isinstance(out_handle, TextIOWrapper):
            out_handle.write(line)

