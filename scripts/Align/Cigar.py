import re
from typing import List

"""Example:
>>> raw_cigar = '4I1=12I3=3I1=2I1=1I1=1I347='
>>> parsed_cigar = Cigar.fromstring(raw_cigar)
>>> parsed_cigar
Cigar(
    cigar_raw='4I1=12I3=3I1=2I1=1I1=1I347=', 
    cigar_str=['I', '=', 'I', '=', 'I', '=', 'I', '=', 'I', '=', 'I', '='], 
    cigar_int=[4, 1, 12, 3, 3, 1, 2, 1, 1, 1, 1, 347]
    )
>>> 
>>> sum(parsed_cigar.cigar_int)
377
>>> parsed_cigar.nth_cigar(-1)
(30, 377)
>>> parsed_cigar.where("I=")
[(0, 4), (4, 5), (5, 17), (17, 20), (20, 23), (23, 24), (24, 26), (26, 27),
(27, 28), (28, 29), (29, 30), (30, 377)]

"""
class Cigar(object):
    def __init__(
            self, 
            cigar_raw:str,
            cigar_str:List[str], 
            cigar_int:List[int]
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


    def __iter__(self):
        for i in range(len(self.cigar_int)):
            yield (self.cigar_str[i], self.cigar_int[i])
        return

    @staticmethod
    def fromstring(raw_cigar:str) -> object:
        p_str = re.compile("[M|I|D|N|S|H|P|=|X]")
        p_int = re.compile(r"\d+")
        m_str = p_str.findall(raw_cigar)
        m_int = [int(i) for i in p_int.findall(raw_cigar)]

        return Cigar(raw_cigar, m_str, m_int)


    # get the start/end locations of the cigar provided by `pattern`
    def positions(self, pattern: str) -> List[tuple]:
        pos = []
        for i, idx in zip(self.cigar_str, range(len(self.cigar_str))):
            if i in pattern:
                pos.append(self.nth_cigar(idx))
        return pos


    # get the start/end location of the nth cigar string
    def nth_cigar(self, n:int) -> tuple([int, int]):
        # covert -1 to len(self.cigar_int)
        if n < 0:
            n = len(self.cigar_int) + n
        # base stage
        if n == 0:
            return (0, self.cigar_int[n])
        else:
            # recursive
            start = self.nth_cigar(n - 1)[-1]
            end = start + self.cigar_int[n]
            return (start, end)

    
    # get index of the cigar provided by `pattern`
    def index(self, pattern:str) -> List[int]:
        index = []
        for i, idx in zip(self.cigar_str, range(len(self.cigar_str))):
            if i in pattern:
                index.append(idx)
        return index
    

    # get number of the cigar provided by `pattern`
    def count(self, pattern:str) -> int:
        counts = 0
        for s, i in zip(self.cigar_str, self.cigar_int):
            if s in pattern:
                counts += i
        return counts


# # Test
# raw_cigar = '4I1=12I3=3I1=2I1=1I1=1I347='
# parsed_cigar = Cigar.from_raw(raw_cigar)
# parsed_cigar

# sum(parsed_cigar.cigar_int)
# parsed_cigar.nth_cigar(-1)
# parsed_cigar.where("I=")


