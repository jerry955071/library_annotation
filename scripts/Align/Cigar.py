import re
from typing import List, Tuple

"""
Example:
# parse `raw_cigar` into Cigar() object
>>> raw_cigar = '20S22=12I33=20D21=30S'
>>> parsed_cigar = Cigar.fromstring(raw_cigar)
>>> parsed_cigar
Cigar(cigar_raw='20S22=12I33=20D21=30S', cigar_str='S=I=D=S', cigar_int=(20, 22,
 12, 33, 20, 21, 30))

# iterate over Cigar() object
>>> for i in parsed_cigar:
...     print(i)
('S', (0, 20))
('=', (20, 42))
('I', (42, 54))
('=', (54, 87))
('D', (87, 87))
('=', (87, 108))
('S', (108, 138))

# subscripting Cigar() object
>>> parsed_cigar.nth_cigar(0)
('S', (0, 20))
>>> parsed_cigar.nth_cigar(-1)
('S', (108, 138))

# get query positions of the requested set of cigar character 
>>> parsed_cigar.positions("=")
[(20, 42), (54, 87), (87, 108)]
"""

class Cigar(object):
    def __init__(
            self, 
            cigar_raw:str,
            cigar_str:str, 
            cigar_int:Tuple[int]
        ) -> None:
        # check length
        if not len(cigar_str) == len(cigar_int):
            raise Exception
        # init object
        super().__init__()
        self.cigar_raw = cigar_raw
        self.cigar_str = cigar_str
        self.cigar_int = cigar_int

    
    def __repr__(self) -> str:
        return "Cigar(cigar_raw=%r, cigar_str=%r, cigar_int=%r)" % (
            self.cigar_raw, self.cigar_str, self.cigar_int
        )

    def __str__(self) -> str:
        return self.cigar_raw


    @staticmethod
    def fromstring(raw_cigar:str) -> object:
        p_str = re.compile("[M|I|D|N|S|H|P|=|X]")
        p_int = re.compile(r"\d+")
        m_str = "".join(p_str.findall(raw_cigar))
        m_int = tuple([int(i) for i in p_int.findall(raw_cigar)])

        return Cigar(raw_cigar, m_str, m_int)


    # get the query start/end locations of the cigar requested by `pattern`
    # def positions(self, pos0: int=0, on_query: bool=True) -> List[tuple]:
    #     pos = []
    #     for i, idx in zip(self.cigar_str, range(len(self.cigar_str))):
    #         pos.append((i, self.pos_nth_cigar(idx, pos0, on_query)))
    #     return pos
    

    #
    # def nth_cigar(self, n:int) -> Tuple[str, Tuple[int, int]]:
    #     return (self.cigar_str[n], self.pos_nth_cigar(n))


    # get the query/reference start, end location (0-based) of the nth cigar string
    def positions(self, pos0:int, on_query:bool) -> tuple([int, int]):
        SKIP = "DH" if on_query else "I"
        ret = []
        for i, s in zip(self.cigar_int, self.cigar_str):
            # base stage
            if len(ret) == 0:
                # start 
                if not on_query:
                    start = pos0 if s not in "HS" else pos0 - i
                else:
                    start = pos0
                # end
                if s in SKIP:
                    end = start
                else:
                    end = start + i
                
                ret.append(
                    (s, (start, end))
                    )
            else:
                # start
                start = ret[-1][-1][-1]
                # end 
                if s in SKIP:
                    end = start
                else:
                    end = start + i
                ret.append(
                    (s, (start, end))
                    )
        return ret

    
    # # get index of the cigar provided by `pattern`
    # def index(self, pattern:str) -> List[int]:
    #     index = []
    #     for i, idx in zip(self.cigar_str, range(len(self.cigar_str))):
    #         if i in pattern:
    #             index.append(idx)
    #     return index
    

    # get number of the cigar provided by `pattern`
    def count(self, pattern:str) -> int:
        counts = 0
        for s, i in zip(self.cigar_str, self.cigar_int):
            if s in pattern:
                counts += i
        return counts

    
    def expand(self) -> str:
        out = ""
        for s, i in zip(self.cigar_str, self.cigar_int):
            out += s * i
        return out


# Test
if __name__ == "__main__":
    raw_cigar = '20S22=12I33=20D21=30S'
    # parse `raw_cigar` into Cigar() object
    parsed_cigar = Cigar.fromstring(raw_cigar)
    parsed_cigar

    # subscripting Cigar() object
    # parsed_cigar.nth_cigar(0)
    # parsed_cigar.nth_cigar(-1)
    parsed_cigar.positions(pos0=0, on_query=False)

