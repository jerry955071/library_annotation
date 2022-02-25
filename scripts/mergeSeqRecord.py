from Bio.SeqRecord import SeqRecord
from typing import List
from collections import deque
import hashlib


# try merge 2 SeqRecord objects
def tryMerge(
    forward: SeqRecord, 
    reverse: SeqRecord, 
    w: int, 
    k: int,
    min_overlap=30
    ) -> None:
    align_col = ovlp(forward.seq, reverse.seq, w, k)
    if align_col[2] == -1:
        return None

    if align_col[2] < min_overlap:
        return None # early return if length of the ovelapped region < min_overlap

    else:
        ref, qry, col = align_col
        # TODO:
        return SeqRecord(
            seq=forward.seq[:ref - col] + \
                reverse.seq[qry - col:],
            id=forward.id,
            name=forward.name,
            description="",
            letter_annotations={
                "phred_quality": qscore(forward)[:ref - col] + \
                    qscore(reverse)[qry - col:]
                }
            )


# get qscore from SeqRecord object
def qscore(record:SeqRecord) -> List[int]:
    return record.letter_annotations["phred_quality"]


# overlapping sequences
def ovlp(seq1:str, seq2:str, w:int, k:int) -> dict:
    """
    Perform fast, approximate overlapping of 
    `seq1` and `seq2` using minimizer methods
    """
    minimizers1 = minimizer.minimizer_sampling(seq1, w, k)
    minimizers2 = minimizer.minimizer_sampling(seq2, w, k)
    return getMaxChain(minimizers1, minimizers2, k)
    

# Generator function yielding indice of the matching minimizers
def getMaxChain(m_ref, m_query, k):
    """
    Input:
        m_ref, m_query: minimizers of seq1 and seq2
        
    Output:
        tuple(x, y, w)
        x: end of seq1
        y: end of seq2
        w: [end - start]
    """
    max_anchor = (-1, -1, -1)
    for idx_q in range(len(m_query)):
        mq = m_query[idx_q]
        for idx_r in range(len(m_ref)):
            mr = m_ref[idx_r]
            if mq[1] == mr[1]:
                qstart = mq[0]
                rstart = mr[0]
                n = 0
                while mq[1] == mr[1]:
                    n += 1
                    try:
                        mq = m_query[idx_q + n]
                        mr = m_ref[idx_r + n]
                    except IndexError:
                        break
                qend = m_query[idx_q + n - 1][0] + k
                rend = m_ref[idx_r + n - 1][0] + k
                w = qend - qstart
                if w > max_anchor[2]:
                    max_anchor = (rend, qend, w)

    return max_anchor


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


# for test
from Bio import SeqIO
in1 = "./output_fastp/forward_passed.fq"
in2 = "./output_fastp/reverse_passed.fq"
record_dict1 = SeqIO.index(in1, "fastq")
record_dict2 = SeqIO.index(in2, "fastq")
fwd = record_dict1["01_A07"]
rev = record_dict2["01_A07"].reverse_complement()
tryMerge(fwd, rev, 13, 8)

