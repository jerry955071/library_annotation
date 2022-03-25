from Align.AlignRecord import AlignRecord
from Align.AlignHeader import AlignHeader
from Align.Cigar import Cigar

class AlignIO:
    header_parser = {
        "@HD": AlignHeader.HDparser,
        "@SQ": AlignHeader.SQparser,
        "@RG": AlignHeader.RGparser,
        "@PG": AlignHeader.PGparser,
        "@CO": AlignHeader.COparser
    }

    def __init__(self, handle) -> object:
        try:
            self.handle = open(handle, "r")
        except TypeError:
            self.handle = handle
        self.header = []
    
    def AlignGenerator(self):
        # skip header lines
        for line in self.handle:
            if line[:3] in AlignIO.header_parser.keys():
                header_tag = line[1:3]
                # use rstrip to remove '\n' at the end of line
                line = line.rstrip()
                # get header parser function
                parser = AlignIO.header_parser[line[:3]]
                self.header.append(parser(line))
            else:
                break

        yield AlignIO.AlignParser(line)

        for line in self.handle:
            yield AlignIO.AlignParser(line)
    
    @staticmethod
    def AlignParser(x):
        typ_funcs = {
                "A": str,
                "f": float,
                "i": int,
                "Z": str
                # "H": <TODO>
                # "B": <TODO>
            }
        # use rstrip to remove '\n' at the end of line
        line = x.rstrip()
        splitted = line.split("\t")
        # record optional fields
        opts_dict = {}
        while len(splitted) > 11:
            pop_out = splitted.pop() # TAG:TYPE:VALUE
            tag = pop_out[:2] # the first two character before the first ":"
            typ = pop_out[3:4] # single character between the ":"
            val = pop_out[5:]
            opts_dict[tag] = typ_funcs[typ](val)

        qname, flag, rname, pos, mapq, cigar, \
            rnext, pnext, tlen, seq, qual = splitted
        
        return AlignRecord(
            qname = qname,
            flag = int(flag),
            rname = rname,
            pos = int(pos),
            mapq = int(mapq),
            cigar = Cigar.fromstring(cigar),
            rnext = rnext,
            pnext = pnext,
            tlen = int(tlen),
            seq = seq,
            qual = qual,
            opts = opts_dict
        )

# for test:
if __name__ == '__main__':
    fin = "/home/woodformation/SSD3/CCC/project_Gilead/library_annotation/output_minimap2/merge.sam"
    with open(fin) as handle:
        for align in AlignIO(handle).AlignGenerator():
            print(align)
