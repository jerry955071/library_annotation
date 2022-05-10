"""Represent a Alignment Record"""
from typing import Tuple
from Align.Cigar import Cigar

class AlignRecord:
    def __init__(
        self,
        qname:str,
        flag:int,
        rname:str,
        pos:int,
        mapq:int,
        cigar:Cigar,
        rnext:str,
        pnext:int,
        tlen:int,
        seq:str,
        qual:str,
        opts:dict=None
    ):
        """Create a AlignRecord.
        Arguments:
        - qname - Query template name (string)
        - flag  - bitwise flag (integer)
        - rname - Reference sequence name (string)
        - pos   - 1-based leftmost mapping position (integer)
        - mapq  - Mapping quality (integer)
        - cigar - Cigar string (Cigar)
        - rnext - Reference name of the mate/next read (string)
        - pnext - Position of the mate/next read (integer)
        - tlen  - Observed template length (integer)
        - seq   - Sequence (string)
        - qual  - Quality scores (string)
        - opts  - Dictionary of key-value pairs
        """
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.opts = opts

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(qname={self.qname!r}, flag={self.flag!r}, "
            f"rname={self.rname!r}, pos={self.pos!r}, mapq={self.mapq!r}, "
            f"cigar={self.cigar!r}, rnext={self.rnext!r}, "
            f"pnext={self.pnext!r}, tlen={self.tlen!r}, seq={self.seq!r}, "
            f"qual={self.qual!r}, opts={self.opts!r}"
        )

    def __str__(self) -> str:
        str_out = ""
        for k, v in vars(self).items():
            if not k == "opts":
                str_out += str(v) + "\t"
            else:
                for tag, val in self.opts.items():
                    if type(val) is str:
                        if " " in val:
                            str_out += "%s:Z:%s\t" % (tag, val)
                        else:
                            str_out += "%s:A:%s\t" % (tag, val)
                    if type(val) is int:
                        str_out += "%s:i:%d\t" % (tag, val)
                    if type(val) is float:
                        str_out += "%s:f:%.1f\t" % (tag, val)

        str_out = str_out.rstrip("\t")
        str_out += "\n"        
        
        return str_out


    # calculate values
    def calculate(self, k):
        methods = {
            "align_bln": AlignRecord.get_align_bln,
            "map_col": AlignRecord.get_map_col,
            "n_error": AlignRecord.get_n_error
        }
        # print all available methods
        if k == "help":
            out_message = "Available methods:"
            for i in methods.keys():
                out_message += " %r" % i
            print(out_message)
            return

        try:
            func = methods[k]
            return func(self)
        except KeyError:
            out_message = "Available methods:"
            for i in methods.keys():
                out_message += " %r" % i
            print(out_message)

    
    # calculate align block length
    def get_align_bln(self) -> int:
        align_blen = 0
        for s, i in zip(self.cigar.cigar_str, self.cigar.cigar_int):
            if s in "=XD":
                align_blen += i
        return align_blen


    # get map column (1-based)
    def get_map_col(self) -> Tuple[int]:
        map_start = self.pos
        map_end = map_start + self.get_align_bln() - 1
        return (map_start, map_end)
    

    # get numbers of error
    def get_n_error(self) -> int:
        return self.cigar.count("IXD")

    
    # generte 'N' substituded sequence
    def n_substitued(self, M=None) -> str:
        if not M:
            M = {
                "S": "",
                "H": "",
                "I": "N",
                "D": "N",
                "X": "N"
            }
        n_subed = ""
        for cig, pos in self.cigar.positions(0, True):
            if cig in M.keys():
                if M[cig] != "":
                    n_subed += M[cig] * (pos[1] - pos[0])
            else:
                n_subed += self.seq[pos[0]:pos[1]]

        return n_subed

    def get_softclip(self) -> tuple:
        front = back = 0
        if self.cigar.cigar_str[0] == "S":
            front = self.cigar.cigar_int[0]
        if self.cigar.cigar_str[-1] == "S":
            back = self.cigar.cigar_int[-1]
        return (front, back)