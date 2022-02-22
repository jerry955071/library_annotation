from Bio.SeqRecord import SeqRecord
from typing import List
from collections import deque
import hashlib


# try merge 2 SeqRecord objects
def tryMerge(forward: SeqRecord, reverse: SeqRecord, min_overlap=30) -> None:
    align_col = ovlp(forward.seq, reverse.seq)
    if align_col["seq1"][0] == -1:
        return None

    align_col_len = align_col["seq1"][1] - align_col["seq1"][0]
    if align_col_len < min_overlap:
        return None # early return if length of the ovelapped region < min_overlap

    else:
        return SeqRecord(
            seq=forward.seq[:align_col["seq1"][0]] + \
                reverse.seq[align_col["seq2"][0]:],
            id=forward.id,
            name=forward.name,
            description="",
            letter_annotations={
                "phred_quality": qscore(forward)[:align_col["seq1"][0]] + \
                    qscore(reverse)[align_col["seq2"][0]:]
                }
            )


# get qscore from SeqRecord object
def qscore(record:SeqRecord) -> List[int]:
    return record.letter_annotations["phred_quality"]


# overlapping sequences
def ovlp(seq1:str, seq2:str, w:int=13, k:int=8) -> dict:
    """
    Perform fast, approximate overlapping of 
    `seq1` and `seq2` using minimizer methods
    """
    minimizers1 = minimizer.minimizer_sampling(seq1, w, k)
    minimizers2 = minimizer.minimizer_sampling(seq2, w, k)
    max_chain = getMaxChain(minimizers1, minimizers2)
    if len(max_chain) == 0:
        return {"seq1": (-1, -1), "seq2": (-1, -1)}
    else:
        return {
            "seq1": (
                minimizers1[max_chain[0][0]][0], 
                minimizers1[max_chain[1][0]][0]
                ),
            "seq2": (
                minimizers2[max_chain[0][1]][0],
                minimizers2[max_chain[1][1]][0]
                )
        }


# get indice of the longest minimizer chain
def getMaxChain(
        m1:List[tuple[int, int]], 
        m2:List[tuple[int, int]]
        ) -> List[tuple[int, int]]:
    """
    Input:
        m1, m2: minimizers of seq1 and seq2
    Output:
        Indices of the longest minimizer chain.
        Output `[]` if no minimizer chain is found.
    """
    chain = []
    max_chain = []
    max_len = 0
    for anchor in anchorGenerator(m1, m2):
        if len(chain) == 0:
            chain.append(anchor)
            len_chain = 1

        elif len(chain) == 1:
            chain.append(anchor)
            if chain[1][0] - chain[0][0] != 1 and chain[1][1] - chain[0][1] != 1:
                chain.pop(0)
            len_chain = 2

        else:
            if chain[1][0] - anchor[0] == -1 and chain[1][1] - anchor[1] == -1:
                chain[1] = anchor
                len_chain += 1
            elif len_chain > max_len:
                max_chain = chain
                chain = [anchor]

    if len(max_chain) == 0:
        if len(chain) == 2:
            max_chain = chain

    return max_chain


# Generator function yielding indice of the matching minimizers
def anchorGenerator(m1, m2):
    """
    Generator function yielding indices of the matching minimizers 
    between `m1` and `m2`
    """
    a1 = a2 = 0
    while (a1 <= len(m1) - 1) & (a2 <= len(m2) - 1):
        while m1[a1][1] != m2[a2][1]:
            if a2 < len(m2) - 1:
                a2 += 1
            elif a2 == len(m2) - 1:
                a2 = 0
                if a1 == len(m1) - 1:
                    return
                a1 += 1
        yield a1, a2
        a1 += 1
    return


class minimizer:
    # sample minimizers
    def minimizer_sampling(seq, w, k, method="standard"):
        """
        Modified from the pseudo code:
        https://academic.oup.com/bioinformatics\
            /article/36/Supplement_1/i111/5870473
        """
        kmers = minimizer.get_kmers(seq, k)
        M = []
        Q = deque()
        for i in range(len(kmers)):
            pos, order = i, minimizer.hash_func(kmers[i])
            while not len(Q) == 0 and Q[-1][1] > order:
                Q.pop() # pop_back()
            Q.append((pos, order))
            if Q[0][0] <= i - w:
                while Q[0][0] <= i - w:
                    Q.popleft() # remove k-mers that slipped out of window
                if method == "winnowing":
                    minimizer.furtherPop(Q) # for winnowing algorithm
            minimizer.sample(Q[0], M)
            if method == "standard":
                minimizer.furtherSample(Q, M) \
                    # for standard minimizer sampling algorithm
        return M


    def sample(minn, M):
        if len(M) == 0:
            M.append(minn)
        else:
            if M[-1] != minn:
                M.append(minn)


    def getOrder(mer):
        return minimizer.hash_func(mer, 8)


    def furtherPop(Q):
        # for winnowing minimizer sampling
        while len(Q) >= 2 and Q[0][1] == Q[1][1]:
            Q = Q[1:]
        return 


    def furtherSample(Q, M):
        # for standard minimizer sampling
        while len(Q) >= 2 and Q[0][1] == Q[1][1]:
            Q.pop()
            minimizer.sample(Q[0], M)
        return


    def get_kmers(seq:str, k) -> List[str]:
        kmers = []
        while len(seq) >= k:
            kmers.append(seq[:k])
            seq = seq[1:]
        return kmers


    def hash_func(string:str, n=8) -> int:
        """
        Return a `n-digit` integer replacement of the input `string`
        """
        return int(hashlib.sha256(string.encode('utf-8')).hexdigest(), 16) \
            % 10**n


# # for test
# from Bio import SeqIO
# in1 = "../output_fastp/forward_passed.fq"
# in2 = "../output_fastp/reverse_passed.fq"
# record_dict1 = SeqIO.index(in1, "fastq")
# record_dict2 = SeqIO.index(in2, "fastq")
# fwd = record_dict1["01_B05"]
# rev = record_dict2["01_B05"].reverse_complement()
# tryMerge(fwd, rev)

