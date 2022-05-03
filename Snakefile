from typing import List

# modify config file before run
configfile: "configs/config.json"


# join list of strings with 'space'
def space_join(x:List[str]) -> str:
    return " ".join(x)


# all
REQUSTED_OUTPUTS = ["merged_fastq", "gapped_fastq"]
rule all:
    input:
        expand("outputs/samtools/{statistics}/{fname}.{statistics}", 
            fname=REQUSTED_OUTPUTS,
            statistics=["stats", "flagstat"]
        ),
        expand([
            "outputs/mappingReporter/{fname}_summary.tsv",
            "outputs/mappingReporter/{fname}_hist.txt",
            "outputs/mappingReporter/{fname}_unmapped.fq"
            ], 
            fname=REQUSTED_OUTPUTS
        )
    params:
        logs_dir="./logs"
    shell:
        """
        # # remove empty logs
        # DIR={params.logs_dir}

        # for file in $(ls $DIR); 
        # do
        #     if [[ ! -s $DIR/$file ]]; then
        #         rm $DIR/$file
        #     fi
        # done
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
        
        
        # combine human and mouse cds.fa
        cat {output.human} {output.mouse} > {output.hnm}
        """


# get fastq from sanger sequecing results 
rule sanger2fastq:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        config["raw_sanger"]
    output:
        directory("outputs/raw_fastq")
    log:
        "logs/sanger2fastq.log"
    shell:
        """
        python scripts/sanger2fastq.py \
            --in_dir {input} \
            --out_dir {output} \
            2> {log} 1> {log}
        """


# merge fastq by sample name 
rule mergeFastq:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        "outputs/raw_fastq/"
    output:
        merged=directory("outputs/merged_fastq/"),
        gapped=directory("outputs/gapped_fastq/")
    params:
        space_join(config["mergeFastq"]["params"])
    log:
        "logs/mergeFastq.log"
    shell:
        """
        python scripts/mergeFastq.py \
            {params} \
            --in_dir {input} \
            --merged_dir {output.merged} \
            --gapped_dir {output.gapped} \
            2> {log} 1> {log}
        """

# indexing reference
rule minimap2_index:
    threads: 8
    input:
        "genomic_data/{target}.fa"
    output:
        "genomic_data/{target}.mmi"
    params:
        space_join(config["minimap2_index"]["params"])
    log:
        "logs/minimap2_index_{target}.log"
    shell:
        """
        minimap2 \
            {params} \
            -d {output} \
            {input} \
            2> {log} 1> {log}
        """


# align reads to human and mouse protein coding transcripts
rule minimap2:
    threads: 8
    input:
        in_dir="outputs/{fastq_dir}",
        ref_idx="genomic_data/gencode.HM.cds.mmi"
    output:
        "outputs/minimap2/{fastq_dir}.sam"
    params:
        space_join(config["minimap2"]["params"])
    log:
        "logs/minimap2_{fastq_dir}.log"
    shell:
        """
        cat {input.in_dir}/*.f*q | # cat all fastq/fq files
        minimap2 \
            {params} \
            -o {output} \
            {input.ref_idx} \
            - \
            2> {log} 1> {log}
        """


# generate stats
rule samtools_stats:
    input:
        "outputs/minimap2/{fname}.sam"
    output:
        "outputs/samtools/stats/{fname}.stats"
    log:
        "logs/samtools_stats_{fname}.log"
    shell:
        """
        samtools stats \
            {input} > {output} \
            2> {log}
        """


# generate flagstat
rule samtools_flagstat:
    input:
        "outputs/minimap2/{fname}.sam"
    output:
        "outputs/samtools/flagstat/{fname}.flagstat"
    log:
        "logs/samtools_flagstat_{fname}.log"
    shell:
        """
        samtools flagstat \
            {input} > {output} \
            2> {log}
        """


# remove secondary/supplementary records
rule samtools_view:
    input:
        "outputs/minimap2/{fname}.sam"
    output:
        "outputs/samtools/filtered/{fname}.filtered.sam"
    params:
        space_join(config["samtools_view"]["params"])
    log:
        "logs/samtools_view_{fname}.log"
    shell:
        """
        samtools view \
            {params} \
            -o {output} \
            {input} \
            2> {log} 1> {log}
        """


# design primer
Primer3 = "tools/primer3_core"
p3_settings_file = "configs/primer3.config"
rule primer3_caller:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        gapped_fq="outputs/gapped_fastq/",
        filtered_sam="outputs/samtools/filtered/gapped_fastq.filtered.sam"
    output:
        directory("outputs/primer3_caller/")
    log:
        stdout="logs/primer3_caller.out",
        stderr="logs/primer3_caller.err"
    shell:
        """
        python scripts/primer3_caller.py \
            --p3 {Primer3} \
            --p3_settings_file {p3_settings_file} \
            --gapped_fq {input.gapped_fq} \
            --filtered_sam {input.filtered_sam} \
            --out_dir {output} \
            2> {log.stderr} 1> {log.stdout}
        """
        

# generate mapping reports
rule mappingReporter:
    conda:
        "envs/project_gilead_chen.yml"
    input:
        "outputs/samtools/filtered/{fname}.filtered.sam"
    output:
        out_summary="outputs/mappingReporter/{fname}_summary.tsv",
        out_hist="outputs/mappingReporter/{fname}_hist.txt",
        out_unmapped="outputs/mappingReporter/{fname}_unmapped.fq"
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
            "outputs/mappingReporter/{fname}_unmapped.fq",
            fname=["merged", "paired1", "paired2", "single"]
            )
    output:
        fa="outputs/blastUnmapped/unmapped.fa",
        blst="outputs/blastUnmapped/blast.out"
    params:
        blastn=space_join(
            config["blastUnmapped"]["params"]["blastn"]
            ),
        out_dir="outputs/blastUnmapped/"
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
            {params.blastn} \
            -query {output.fa} \
            -out {output.blst} \
            2> {log.blst} 1> {log.blst}

        # split blast results
        python scripts/splitBlast.py \
            -i {output.fa} \
            --out_dir {params.out_dir} \
            2> {log.py} 1> {log.py}
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
        echo {rule}
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


# SKIP THIS STEP
# # fastp was use to perform:
#     # 1. quality trimming
#     # 2. convert phred 64 to phred 33
#     # 3. generate an overall report    
# rule fastp_qc:
#     input:
#         "raw_fastq/{fname}.fq"
#     output:
#         out="output_fastp/{fname}_passed.fq",
#         failed="output_fastp/{fname}_failed.fq",
#         report_json="output_fastp/{fname}_report.json",
#         report_html="output_fastp/{fname}_report.html"
#     params:
#         space_join(config["fastp_qc"]["params"])
#     log:
#         "logs/fastp_{fname}.log"
#     shell:
#         """
#         fastp \
#             {params} \
#             --in1 {input} \
#             --out1 {output.out} \
#             --failed_out {output.failed} \
#             --json '{output.report_json}' \
#             --html '{output.report_html}' \
#             2> {log} 1> {log}
#         """
