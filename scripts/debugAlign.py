from Align.AlignIO import AlignIO

file = AlignIO("../output_minimap2/merged_human.sam").AlignGenerator()
for line in file:
    pass