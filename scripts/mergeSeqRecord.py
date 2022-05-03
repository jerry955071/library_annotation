from Bio.SeqRecord import SeqRecord
from Align.Cigar import Cigar
from typing import List
from collections import deque, namedtuple
import hashlib
import edlib


# NamedTuple recording between-sequence information
Anchor = namedtuple("Anchor", "rstart, rend, qstart, qend")
Chain = namedtuple("Chain", "rstart, rend, qstart, qend")

# try merge 2 SeqRecord objects
def tryMerge(
    forward: SeqRecord, 
    reverse: SeqRecord, 
    w: int, 
    k: int,
    min_overlap=30
    ) -> None:
    
    # get maximized minimizer chain
    max_chain = getMaxChain(ref=forward.seq, qry=reverse.seq, w=w, k=k)
    
    # return `None` if no chain found
    if any([i == -1 for i in max_chain]):
        return None

    # return `None` if ovelapped length < min_overlap
    if max_chain.rend - max_chain.rstart < min_overlap:
        return None

    else:
        # perform global pairwise alignment at chained region
        sub_qry = reverse[max_chain.qstart:max_chain.qend]
        sub_ref = forward[max_chain.rstart:max_chain.rend]
        edlib_align = edlib.align(
            query=str(sub_qry.seq),
            target=str(sub_ref.seq),
            mode="NW",
            task="path"
        )

        # for each column in the alignment
        generator_qq = (i for i in qscore(sub_qry))
        generator_qr = (i for i in qscore(sub_ref))
        generator_bq = (i for i in sub_qry)
        generator_br = (i for i in sub_ref)
        combined_bases = ""
        combined_quals = []
        last_qq = deque(maxlen=5) # record up 5 previous q-score
        last_qr = deque(maxlen=5) # record up 5 previous q-score
        for cig in Cigar.fromstring(edlib_align["cigar"]).expand():
            bq = next(generator_bq) if cig != "D" else "-"
            qq = next(generator_qq) if cig != "D" else round(sum(last_qq)/len(last_qq))
            br = next(generator_br) if cig != "I" else "-"
            qr = next(generator_qr) if cig != "I" else round(sum(last_qr)/len(last_qr))
            last_qq.append(qq)
            last_qr.append(qr)
            
            # pick the higher-q base 
            base = br if qr > qq else bq
            # pick the higher q
            qual = max(qr, qq)
            # skip "-" characters
            if base == "-":
                continue
            else:
                combined_bases += base
                combined_quals.append(qual)

        return SeqRecord(
            seq=forward.seq[:max_chain.rstart] + \
                combined_bases + \
                reverse.seq[max_chain.qend:],
            id=f"{forward.id}|{reverse.id}",
            name=forward.name,
            description="",
            letter_annotations={
                "phred_quality": \
                    qscore(forward)[:max_chain.rstart] + \
                    combined_quals + \
                    qscore(reverse)[max_chain.qend:]
                }
        )


# get qscore from SeqRecord object
def qscore(record:SeqRecord) -> List[int]:
    return record.letter_annotations["phred_quality"]


def getMaxChain(ref:str, qry:str, w:int, k:int) -> tuple:
    """
    Input:
        m_ref, m_qry: minimizers of seq1 and seq2
        
    Output:
        Chain(rstart, rend, qstart, qend)
    """
    m_ref = minimizer.minimizer_sampling(ref, w, k)
    m_qry = minimizer.minimizer_sampling(qry, w, k)
    max_chain = Chain(-1, -1, -1, -1)

    for idx_q in range(len(m_qry)):
        mq = m_qry[idx_q]
        for idx_r in range(len(m_ref)):
            mr = m_ref[idx_r]
            # found seed
            if mq[1] == mr[1]:
                qstart = mq[0]
                rstart = mr[0]
                n = 0
                # extend
                while mq[1] == mr[1]:
                    n += 1
                    try:
                        mq = m_qry[idx_q + n]
                        mr = m_ref[idx_r + n]
                    except IndexError:
                        break
                qend = m_qry[idx_q + n - 1][0] + k
                rend = m_ref[idx_r + n - 1][0] + k
                
                # evaluate new chain length
                min_len = min(qend - qstart, rend - rstart)
                if min_len > min(
                    max_chain.rend - max_chain.rstart,
                    max_chain.qend - max_chain.qstart
                    ):
                    max_chain = Chain(rstart, rend, qstart, qend)
    
    return max_chain


# un-used
def get_anchors(
    reference_minimizer,
    query_minimizer,
    kmer_size
    ) -> List[tuple]:
    """
    Find list of anchors
    A 'anchor' is defined as tuple(rend, qend, anchor_length), while the
    'anchor_length' is equal to (rend - rstart) and (qend - qstart)
    """
    anchors_list = []
    for m_ref in reference_minimizer:
        for m_qry in query_minimizer: 
            if m_ref[1] == m_qry[1]:
                anchors_list.append(
                    Anchor(
                        rstart=m_ref[0],
                        rend = m_ref[0] + kmer_size,
                        qstart=m_qry[0],
                        qend = m_qry[0] + kmer_size
                        )
                    )
    
    return anchors_list


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
if __name__ == '__main__':
    from Bio import SeqIO
    in1 = "../raw_fastq/forward.fq"
    in2 = "../raw_fastq/reverse.fq"
    record_dict1 = SeqIO.index(in1, "fastq")
    record_dict2 = SeqIO.index(in2, "fastq")
    for i in list(record_dict1.keys()):
        try:
            fwd = record_dict1[i]
            rev = record_dict2[i].reverse_complement()
            tryMerge(fwd, rev, 13, 8)
        except KeyError:
            pass

