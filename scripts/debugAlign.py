from Align.AlignIO import AlignIO
import sys

fin = sys.argv[1]
file = AlignIO(fin).AlignGenerator()
for line in file:
    if line.rname == "*":
        pass