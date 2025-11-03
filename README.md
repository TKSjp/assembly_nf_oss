# assembly_nf

Nextflow assembly pipelines for Entamoeba.

## environment

- This script is written for super computer in SGE system, especially shirokane (https://gc.hgc.jp)
- You need nextflow and conda environment written in config files.

## setup of conda environment

see config and prepare all environments.

## run

here, you need log dir, input dir, and script dir which has each pipeline script.

- qsub -N assembly script/assembly.sh

## output of the pipeline

results are in result directory

- hifiasm
    - (readname).bp.p_ctg.fna is the final contigs.
- minimap2_remap
- quast
- samtools_sort
- trimmed_reads

## contact

If you have any troubles, please feel free to open an issue.
Tetsuro Kawano-Sugaya (sugaya@tetsu.ro)
https://tetsu.ro