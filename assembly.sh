#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -pe def_slot 8
#$ -l s_vmem=8G,mem_req=8G
#$ -q '!mjobs_rerun.q'

# usage: qsub -N assembly script/assembly.sh

source ~/.bashrc.intr
module load java/11
unset JAVA_TOOL_OPTIONS
export JAVA_TOOL_OPTIONS="-XX:+UseG1GC -Xmx6g -Xms6g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=2"
export NXF_VER=23.10.1
export NXF_EXECUTOR=sge
export NXF_OPTS="-Xms6g -Xmx6g"

nextflow -log log/assembly/$(date +%Y%m%d_%H%M%S).txt run script/assembly.nf -c script/assembly.config -profile standard -with-report report/assembly/$(date +%Y%m%d_%H%M%S).html -with-trace trace/assembly/$(date +%Y%m%d_%H%M%S).txt -work-dir work/assembly/$(date +%Y%m%d_%H%M%S) -resume
