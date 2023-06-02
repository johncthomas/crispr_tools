# crispr_tools
Collection of tools for CRISPR-Cas9 genetic screen analysis.

## Introduction
These scripts neatly pipeline read counting from FASTQ and wrap differential analysis with DrugZ and MAGeCK in batch.

## Pipeline overview
CRISPR screen analysis using crispr_tools proceeds by the following steps:
1. From pre-aligned FASTQ, determine read counts per protospacer and map them to the reference CRISPR library.
2. Specify replicate-to-condition mappings in the experimental design and select differential methods to run.
3. Run the analysis in batch.

# Tutorial
A tutorial dataset is provided in the `tut` directory. To verify that you have correctly run the analysis, your directory should resemble `tut_ran`. Follow the steps below for a detailed explanation of the analysis.

## Reading counts
The `count_reads.py` script handles read counting from pre-aligned FASTQ files and mapping to a reference library, if supplied.

`count_reads.py [FASTQ_FILES] -s SLICE --merge-samples --just-go -p OUTDIR --library=LIB`

The slice (`SLICE`) is range that specifies the position of the protospacer along the reads. It is supplied as zero-based half-open indices. FASTQ files must be pre-aligned and have a constant slice across all reads in all runs.

The slice can be determined by inspecting nucleotide frequencies as a function of cycle. It appears as a region of uniform probability across the four bases. The `ntByCycle.R` script in the exorcise package (https://github.com/SimonLammmm/exorcise/) can generate such nucleotide traces.

`ntByCycle.R -f [FASTQ_FILE] -o OUTDIR`

Specify `--merge-samples` if you have technical replicates. You'll know this is the case if you have FASTQ files with `_L001` and `_L002` or similar.

Specify the path of the reference library (`LIB`). It must have `seq`, `gene`, and `guide` columns, indicating protospacer sequence, gene symbol, and unique guide IDs, respectively.

Running the `count_reads.py` script generates in the `OUTDIR` raw counts per protospacer for each FASTQ. It will also generate counts per gene symbol if specified with a reference library.

## Design matrix
Populate the Excel workbook with experimental details pertaining to the run.

In the “Experimental details” sheet, you must supply the analysis version and file prefix. Results will be saved at `file_prefix`/`analysis_version/`. All other fields are for your own convenience should you need to look at this analysis again in the future.

In the “Sample details” sheet, ensure that entries in the “Replicate” column correspond to column names in the counts file. Biological replicates should be given the same value in the “Sample” column. All other columns are optional but are here for your convenience.

In the “Control groups” sheet, ensure that control and test samples defined here match those in the “Sample” column of the “Sample details” sheet. In the “Group” column, use names to group different comparisons together. You can define different settings for each group in the “Analyses” sheet.

In the “Analyses” sheet, ensure that the “Control group” entries match those in the “Group” column of the “Control groups” sheet. For “Paired”, set TRUE if the replicates for the control and test samples with the same index are meaningfully paired; otherwise set FALSE. In “Method”, use drugz, mageck, or drugz,mageck according to which analyses you want to conduct. In “Counts file”, supply the absolute or relative path of the counts file. In “Add pseudocount”, use 5 unless you have a reason to use a different number. The other fields can be left blank.

## Run batch differential analysis
This step is the last step and uses `crispr_pipeline.py`.

`crispr_pipeline.py DET`

Specify as a positional argument the experimental details Excel workbook. You will find in `file_prefix`/`analysis_version` the standard results files from a DrugZ and/or MAGeCK analysis, depending on your selection.

