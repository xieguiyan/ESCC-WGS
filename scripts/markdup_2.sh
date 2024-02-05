#!/bin/bash
sample="${snakemake_wildcards[sample]}"
inputbam="${snakemake_input[SortedBam]}"
if [ $(${snakemake_params[samtools]} view -h ${snakemake_input[SortedBam]} | \
   head -n 20000 |grep "^@RG.*SM:${sample}" |wc -l) -lt 1 ]
then
   inputbam=$(echo "${inputbam}" |sed "s/\.bam$/\.pic\.bam/")
   time ${snakemake_params[java]} -jar ${snakemake_params[picard]} \
      AddOrReplaceReadGroups -I ${snakemake_input[SortedBam]} -O ${inputbam} \
      --RGID other_${sample} --RGLB lib1${sample} \
      --RGPL illumina \
      --RGSM ${sample} --VALIDATION_STRINGENCY LENIENT
fi
prefix=$(echo "${snakemake_output[SortedMarkdupBam]}" |sed "s/\.bam$//")
outputp=$(echo "${prefix}" |grep -o "/.*/")
if [ ! -d ${outputp} ]; then mkdir -p ${outputp} ;fi
time ${snakemake_params[sentieon]} driver -t 20 -i \
   ${inputbam} --algo LocusCollector \
   --fun score_info ${prefix}.score.txt
time ${snakemake_params[sentieon]} driver -t 20 -i \
   ${inputbam} --algo Dedup \
   --score_info ${prefix}.score.txt \
   --metrics ${prefix}.dedup.metrics \
   ${snakemake_output[SortedMarkdupBam]}

[ ! -f ${snakemake_output[SortedMarkdupBam]}.bai ] && \
   ${snakemake_params[samtools]} index ${snakemake_output[SortedMarkdupBam]}
[ ! -f ${snakemake_output[SortedMarkdupBai]} ] && \
   ln -s ${snakemake_output[SortedMarkdupBam]}.bai  ${snakemake_output[SortedMarkdupBai]}
echo "[## markdup ] date: $(date), parameter: ${sample}, status: finished ..." >>${snakemake_config[status]}

