# modify config file before run
configfile: "config.json"


# all
rule all:
    input:
        expand("output_samtools/stats/{fname}.stats", 
            fname=["merged", "single", "paired1", "paired2"],
        ),
        expand("output_samtools/flagstats/{fname}.flagstat", 
            fname=["merged", "single", "paired1", "paired2"]
        ),
        expand("output_mappingReporter/{fname}_summary.tsv", 
            fname=["merged", "single", "paired1", "paired2"]
        ),
        "output_blastUnmapped/blast.out"
    params:
        logs_dir="./logs"
    shell:
        """
        # remove empty logs
        DIR={params.logs_dir}

        for file in $(ls $DIR); 
        do
            if [[ ! -s $DIR/$file ]]; then
                rm $DIR/$file
            fi
        done
        """


# extract CDS from Gencode reference
rule extract_cds:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        human="genomic_data/gencode.v38.pc_transcripts.fa",
        mouse="genomic_data/gencode.vM27.pc_transcripts.fa"
    output:
        human="genomic_data/gencode.v38.cds.fa",
        mouse="genomic_data/gencode.vM27.cds.fa",
        hnm="genomic_data/gencode.HM.cds.fa"
    log:
        "logs/extract_cds.log"
    shell:
        """
        # human pc_transcripts to cds
        python scripts/extract_cds.py \
            -i {input.human} \
            -o {output.human} \
            2> {log} 1> {log}
        
        # mouse pc_transcripts to cds
        python scripts/extract_cds.py \
            -i {input.mouse} \
            -o {output.mouse} \
            2> {log} 1> {log}
        
        
        # combine human and mouse cds
        cat {output.human} {output.mouse} > {output.hnm}
        """
 

# get fastq from sanger sequecing results 
rule sanger2fastq:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        config["input_folder"]
    output:
        out1="raw_fastq/forward.fq",
        out2="raw_fastq/reverse.fq"
    log:
        "logs/sanger2fastq.log"
    shell:
        """
        python scripts/sanger2fastq.py \
            -i {input} \
            --out1 {output.out1} \
            --out2 {output.out2} \
            --suff1 ADF.ab1 \
            --suff2 ADR.ab1 \
            2> {log} 1> {log}
        """


# fastp was use to perform:
    # 1. quality trimming
    # 2. convert phred 64 to phred 33
    # 3. generate an overall report    
rule fastp_qc:
    input:
        "raw_fastq/{fname}.fq"
    output:
        out="output_fastp/{fname}_passed.fq",
        failed="output_fastp/{fname}_failed.fq"
    params:
        length_required=50,
	    mean_q=20,
        cut_window_size=10,
        out_dir="output_fastp"
    log:
        "logs/fastp_{fname}.log"
    shell:
        """
        fastp \
            --phred64 \
            --in1 {input} \
            --out1 {output.out} \
            --failed_out {output.failed} \
            --cut_mean_quality {params.mean_q} \
            --cut_front \
            --cut_tail  \
            --cut_window_size {params.cut_window_size} \
	        --length_required {params.length_required} \
            --json '{params.out_dir}/{wildcards.fname}_report.json' \
            --html '{params.out_dir}/{wildcards.fname}_report.html' \
            2> {log} 1> {log}
        """


# classify QC-passed reads into merged/single/paired1/paired2.fq
rule classifier:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        in1="output_fastp/forward_passed.fq",
        in2="output_fastp/reverse_passed.fq"
    output:
        merged="output_classifier/merged.fq",
        single="output_classifier/single.fq",
        paired1="output_classifier/paired1.fq",
        paired2="output_classifier/paired2.fq",
        report="output_classifier/report.txt",
        summary="output_classifier/summary.tsv"
    params:
        w="13",
        k="8",
        min_ovlp="30"
    log:
        "logs/classifier.log"
    shell:
        """
        python scripts/classifier.py \
            -w {params.w} -k {params.k} --min_ovlp {params.min_ovlp} \
            --in1 {input.in1} \
            --in2 {input.in2} \
            --out_merged {output.merged} \
            --out_single {output.single} \
            --out_paired1 {output.paired1} \
            --out_paired2 {output.paired2} \
            --summary {output.summary} \
            --report {output.report} \
            2> {log} 1> {log}
        """


# align reads to human and mouse protein coding transcripts
rule minimap2:
    threads: 8
    input:
        query="output_classifier/{fname}.fq",
        ref="genomic_data/gencode.HM.cds.fa"
    output:
        "output_minimap2/{fname}.sam"
    log:
        "logs/minimap2_{fname}.log"
    shell:
        """
        minimap2 -ax asm5 \
            -E4,3 `# gap extension penalty` \
            --eqx \
            --no-end-flt `# dont filter end seed before base-level alignment` \
            -t {threads} \
            -o {output} \
            {input.ref} \
            {input.query} \
            2> {log} 1> {log}
        """


# generate stats
rule samtools_stats:
    input:
        "output_minimap2/{fname}.sam"
    output:
        "output_samtools/stats/{fname}.stats"
    log:
        "logs/samtools_stats_{fname}.log"
    shell:
        """
        samtools stats \
            {input} > {output} \
            2> {log}
        """


# generate flag stats
rule samtools_flagstats:
    input:
        "output_minimap2/{fname}.sam"
    output:
        "output_samtools/flagstats/{fname}.flagstat"
    log:
        "logs/samtools_flagstats_{fname}.log"
    shell:
        """
        samtools flagstat \
            {input} > {output} \
            2> {log}
        """


# remove secondary/supplementary records
rule samtools_view:
    input:
        "output_minimap2/{fname}.sam"
    output:
        "output_samtools/filtered/{fname}.sam"
    log:
        "logs/samtools_view_{fname}.log"
    shell:
        """
        samtools view \
            -F 256 `# remove secondary alignments` \
            -F 2048 `# remove supplementary alignments` \
            -h \
            -o {output} \
            {input} \
            2> {log} 1> {log}
        """


# generate mapping reports
rule mappingReporter:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        "output_samtools/filtered/{fname}.sam"
    output:
        out_summary="output_mappingReporter/{fname}_summary.tsv",
        out_hist="output_mappingReporter/{fname}_hist.txt",
        out_unmapped="output_mappingReporter/{fname}_unmapped.fq"
    log:
        "logs/mappingReporter_{fname}.log"
    shell:
        """
        python scripts/mappingReporter.py \
            -i {input} \
            --out_summary {output.out_summary} \
            --out_hist {output.out_hist} \
            --out_unmapped {output.out_unmapped} \
            2> {log} 1> {log}
        """


# collect, blast all unmapped reads
rule blastUnmapped:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        expand(
            "output_mappingReporter/{fname}_unmapped.fq",
            fname=["merged", "paired1", "paired2", "single"]
            )
    output:
        fa="output_blastUnmapped/unmapped.fa",
        blst="output_blastUnmapped/blast.out"
    params:
        max_target_seqs=5,
        db="nt",
        outfmt=0
        out_dir="output_blastUnmapped/"
    log:
        py="logs/collectUnmapped.log",
        blst="logs/blastUnmapped.log"
    shell:
        """
        # collect unmapped reads
        python scripts/collectUnmapped.py \
            -i {input} \
            --out_fa {output.fa} \
            2> {log.py} 1> {log.py}

        # blast unmapped reads
        blastn \
            -query {output.fa} \
            -db {params.db} \
            -remote \
            -max_target_seqs {params.max_target_seqs} \
            -outfmt {params.outfmt} \
            -out {output.blst} \
            2> {log.blst} 1> {log.blst}

        # format blast results
        python scripts/splitBlast.py \
            -i {output.fa} \
            --out_dir {params.out_dir}
        """


# ==============================================================================
# un-used rules:
rule template:
    input:
    output:
    log:
    shell:
        """
        echo Hello world!
        """


# fastp: merge oberlapping paired-end reads
# rule fastp_classify:
#     input:
#         in1="output_fastp/forward_passed.fq",
#         in2="output_fastp/reverse_passed.fq"
#     output:
#         merged="output_fastp_classify/merged.fq",
#         paired1="output_fastp_classify/paired1.fq",
#         paired2="output_fastp_classify/paired2.fq",
#         single="output_fastp_classify/single.fq"
#     params:
#         min_ovlp=30,
#         out_dir="output_fastp_classify"
#     log:
#         "logs/fastp_classify.log"
#     shell:
#         """
#         fastp \
#             --in1 {input.in1} \
#             --in2 {input.in2} \
#             --out1 {output.paired1} \
#             --out2 {output.paired2} \
#             --unpaired1 {output.single} \
#             --unpaired2 {output.single} \
#             --merge \
#             --merged_out {output.merged} \
#             --overlap_len_require {params.min_ovlp} \
#             --json '{params.out_dir}/report.json' \
#             --html '{params.out_dir}/report.html' \
#             2> {log} 1> {log}
#         """


# # align reads to human and mouse protein coding transcripts
# rule minimap2:
#     threads: 4
#     input:
#         query="output_classifier/{fname}.fq",
#         ref_h="genomic_data/gencode.v38.cds.fa",
#         ref_m="genomic_data/gencode.vM27.cds.fa"
#     output:
#         human="output_minimap2/{fname}_human.sam",
#         mouse="output_minimap2/{fname}_mouse.sam"
#     log:
#         human="logs/minimap2_{fname}_human.log",
#         mouse="logs/minimap2_{fname}_mouse.log"
#     shell:
#         """
#         # map to human protein coding transcripts
#         minimap2 -ax asm5 --eqx \
#             -t {threads} \
#             -o {output.human} \
#             {input.ref_h} \
#             {input.query} \
#             2> {log.human} 1> {log.human}

#         # map to mouse protein coding transcripts
#         minimap2 -ax asm5 --eqx \
#             -t {threads} \
#             -o {output.mouse} \
#             {input.ref_m} \
#             {input.query} \
#             2> {log.mouse} 1> {log.mouse}
#         """

