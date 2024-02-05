#!/bin/bash

sample="${snakemake_wildcards[sample]}"
GATK="${snakemake_params[GATK]}"
java="${snakemake_params[java]}"
reference="${snakemake_params[reference]}"
outputp=$(echo "${snakemake_output[completed]}" |grep -o "/.*/")
if [ ! -d ${outputp} ];then mkdir -p ${outputp} ;fi
echo "ignore this part..." >${snakemake_output[ignore]}
