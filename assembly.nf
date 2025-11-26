#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process trim_hifi_reads {
    publishDir "result/trimmed_reads", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("${id}_trim.fastq.gz"), emit: trimmed_fastq
    path "version_${id}.txt"

    shell:
    '''
    seqkit seq -m 5000 -j !{task.cpus} !{fastq} | pigz -p !{task.cpus} > !{id}_trim.fastq.gz
    seqkit -h | head -n 3 > version_!{id}.txt
    '''
}

process minimap2_filter {
    publishDir "result/minimap2_filter", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fastq)
    path(fasta)

    output:
    tuple val(id), path("filter_${id}.sam"), emit: sam
    path "version_${id}.txt"

    shell:
    '''
    minimap2 -t !{task.cpus} -ax map-hifi !{fasta} !{fastq} > filter_!{id}.sam
    minimap2 --version > version_!{id}.txt
    '''
}

process samtools_split {
    publishDir "result/samtools_split", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path("${id}_unmapped.fastq.gz"), emit: unmapped_fastq
    tuple val(id), path("${id}_mapped.fastq.gz"), emit: mapped_fastq
    path "stats_${id}.tsv"
    path "version_${id}.txt"

    shell:
    '''
    samtools view -@ !{task.cpus} -b -f 4 !{sam} | samtools fastq -@ !{task.cpus} - | pigz -p !{task.cpus} > !{id}_unmapped.fastq.gz
    samtools view -@ !{task.cpus} -b -F 4 !{sam} | samtools fastq -@ !{task.cpus} - | pigz -p !{task.cpus} > !{id}_mapped.fastq.gz
    seqkit stat *.fastq.gz > stats_!{id}.tsv
    samtools --version | head -n 1 > version_!{id}.txt
    '''
}

process hifiasm {
    publishDir "result/hifiasm", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("${id}.bp.p_ctg.fna"), emit: fna
    path "${id}"

    shell:
    '''
    mkdir -p !{id}
    hifiasm -f0 -z20 -l0 -o !{id} -t !{task.cpus} !{fastq} 2> log_!{id}.txt
    awk '/^S/{print ">"$2;print $3}' !{id}.bp.p_ctg.gfa > !{id}.bp.p_ctg.fna
    cp *.fna *.gfa *.txt !{id}
    hifiasm -v > !{id}/version_!{id}.txt
    '''
}

process quast {
    publishDir "result/quast", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("quast_${id}")
    path "version_${id}.txt"

    shell:
    '''
    quast !{fasta} -o quast_!{id}
    quast -v > version_!{id}.txt
    '''
}

process minimap2_remap {
    publishDir "result/minimap2_remap", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fastq)
    tuple val(id2), path(fna)

    output:
    tuple val(id), path("${id}_mapped.sam"), emit: sam
    path "version_${id}.txt"

    shell:
    '''
    minimap2 -t !{task.cpus} -ax map-hifi !{fna} !{fastq} > !{id}_mapped.sam
    minimap2 --version > version_!{id}.txt
    '''
}

process samtools_sort {
    publishDir "result/samtools_sort", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path("${id}_mapped.sort.bam")
    path "${id}_mapped.sort.bam.bai"

    shell:
    '''
    samtools view -@ !{task.cpus} -bS !{sam} > !{id}_mapped.bam
    samtools sort -@ !{task.cpus} -o !{id}_mapped.sort.bam !{id}_mapped.bam
    samtools index !{id}_mapped.sort.bam
    '''
}

workflow {
    fastq_files = channel.fromFilePairs('input/*.fastq.gz', size: 1)
    trimmed_fastq_files = trim_hifi_reads(fastq_files)
    hifiasm_assembled = hifiasm(trimmed_fastq_files.trimmed_fastq)
    quast_analysed = quast(hifiasm_assembled.fna)
    minimap2_remapped = minimap2_remap(trimmed_fastq_files.trimmed_fastq, hifiasm_assembled.fna)
    samtools_sorted = samtools_sort(minimap2_remapped.sam)
}
