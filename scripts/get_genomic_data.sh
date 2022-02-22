#!/bin/bash
function wget_and_unzip() {
    #
    DIR=$1
    URL=$2
    BASENAME=${URL##*/}
    LOG1=./logs/download_${BASENAME}.log
    LOG2=./logs/decompress_${BASENAME}.log
    # 
    wget --output-file=$LOG1 --directory-prefix=./$DIR/ $URL && \
        gzip -d ./$DIR/$BASENAME 2>&1 $LOG2 &
}

# Get Human protein coding trasncript sequences --------------------------------
# Source: Gencode
# Species: Human
# Version: Release 38 (GRCh38.p13)

## Content: Protein-coding transcript sequences
## Region: CHR
## Description: 
## Nucleotide sequences of coding transcripts on the reference chromosomes
## Transcript biotypes: protein_coding, nonsense_mediated_decay, non_stop_decay,
## IG_*_gene, TR_*_gene, polymorphic_pseudogene
wget_and_unzip "./genomic_data" \
	"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/\
release_38/gencode.v38.pc_transcripts.fa.gz"

# Get Mouse protein coding trasncript sequences --------------------------------
# Source: Gencode
# Species: Mouse
# Version: Release M27 (GRCm39)

## Content: Protein coding trasncript sequences 
## Region: CHR
## Description: 
## Nucleotide sequences of coding transcripts on the reference chromosomes
## Transcript biotypes: protein_coding, nonsense_mediated_decay, non_stop_decay,
## IG_*_gene, TR_*_gene, polymorphic_pseudogene
wget_and_unzip "./genomic_data" \
	"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/\
release_M27/gencode.vM27.pc_transcripts.fa.gz"


