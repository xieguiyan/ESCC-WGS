#!/bin/bash
sample="${snakemake_wildcards[sample]}"
prefix=$(echo "${snakemake_input[0]}" |grep -o "^/.*/")/$sample
outputp=$(echo "${snakemake_output[BAM]}" |grep -o "/.*/")
input_data_format="${snakemake_config[Inputfm]}"
if [ ! -d ${outputp} ]; then  mkdir -p ${outputp} ;fi
if [[ ${input_data_format} == "BAM" ]] && [ -f ${prefix}.bam ]
then
   mv ${prefix}.bam ${snakemake_output[BAM]}
elif [ $(ls ${prefix}.clean*fastq.gz 2>/dev/null |wc -l) -gt 0 ]
then
   metasfastq=""
   if [ -f ${prefix}.clean.1.fastq.gz -a -f ${prefix}.clean.2.fastq.gz ]
   then
      metasfastq="${prefix}.clean.1.fastq.gz ${prefix}.clean.2.fastq.gz"
   elif [ -f ${prefix}.clean.fastq.gz ]
   then
      metasfastq="${prefix}.clean.fastq.gz"
   fi
   Mparam=""
   if ${snakemake_params[marksecondary]}; then echo "bwa \-M carry on"; Mparam="-M" ;fi
   Kparam=""
   if [ ${snakemake_params[Kvalue]} -gt 0 ];then Kparam="-K ${snakemake_params[Kvalue]}" ;fi
   echo "mksy: ${snakemake_params[marksecondary]}, ${snakemake_params[Kvalue]}"
   ${snakemake_params[bwa]} mem -t ${snakemake[threads]} ${Mparam} \
       ${Kparam} \
       -R "@RG\tID:foo_lane_${sample}\tPL:illumina\tLB:library_${sample}\tSM:${sample}" \
       ${snakemake_params[refer]} ${metasfastq} | \
       ${snakemake_params[samtools]} view -Sh -b -@ ${snakemake[threads]} - > \
       ${snakemake_output[BAM]}
else
   echo "ERROR: data not found for bwa... "
   exit -1
fi
time ${snakemake_params[samtools]} sort -@ 4 -m 4G \
   -O bam -o ${snakemake_output[SortedBam]} ${snakemake_output[BAM]}
${snakemake_params[samtools]} index ${snakemake_output[SortedBam]}
echo "[## bwa ] date: $(date), parameter: ${sample}, status: finished ..." >${snakemake_config[status]}
