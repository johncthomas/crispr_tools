# Find the slice in the FASTQ file (requires https://github.com/SimonLammmm/exorcise)
ntByCycle.R -f fastq/Givpir_3.fastq.gz -o ntByCycle/

# Count reads without alignment
mkdir cts
count_reads.py fastq/*.fastq.gz -s 6,24 --just-go -p cts/cts --library=lib/library.tsv

# Perform batch differential analysis
crispr_pipeline.py det/design.xlsx