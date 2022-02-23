from Align.AlignIO import AlignIO

file = AlignIO("../output_samtools/filtered/merged.sam").AlignGenerator()
for line in file:
    if line.rname == "*":
        pass