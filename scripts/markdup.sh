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

time ${snakemake_params[java]} -jar ${snakemake_params[picard]} \
   MarkDuplicates REMOVE_DUPLICATES=true \
   ASSUME_SORTED=true \
   CREATE_INDEX=true \
   I=${inputbam} \
   O=${snakemake_output[SortedMarkdupBam]} \
   M=$(echo ${snakemake_output[SortedMarkdupBam]}| sed "s/\.bam$//")_metrics.txt


[ ! -f  ${snakemake_output[SortedMarkdupBai]} ] && \
   ln -s ${snakemake_output[SortedMarkdupBam]}.bai  ${snakemake_output[SortedMarkdupBai]}
echo "[## markdup ] date: $(date), parameter: ${sample}, status: finished ..." >>${snakemake_config[status]}

