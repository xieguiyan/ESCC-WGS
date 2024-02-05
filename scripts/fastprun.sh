#!/bin/bash
sample="${snakemake_wildcards[sample]}"
outprefix=$(echo ${snakemake_output[0]} |sed "s/[^/]*$//")/${sample}
if [ ! -d ${outputfastq} ] ;then mkdir -p  ${outputfastq} ;fi
wkrawdata=$(echo "${snakemake_input[0]}" | grep -o "/.*/")
inprefix=${wkrawdata}/${sample}
qtprefix=$(echo "${snakemake_input[completed2]}" |grep -o "/.*/" )/${sample}

if [ ! -f ${qtprefix}.fastqc.PASS ]
then
   fastqpara=""
   if [[ $(ls ${wkrawdata} |grep "${sample}_[12]\.fastq\.gz$" |wc -l) -eq 2 ]]
   then
      fastqpara="-i ${inprefix}.1.fastq.gz -o ${outprefix}.clean.1.fastq.gz -I ${inprefix}.2.fastq.gz -O ${outprefix}.clean.2.fastq.gz --detect_adapter_for_pe"
   elif [ -f ${inprefix}.fastq.gz ]
   then
      fastqpara="-i ${inprefix}.fastq.gz -o ${outprefix}.clean.fastq.gz"
   fi
   ${snakemake_params[fastp]} -w 12 $fastqpara \
      -q 15 -u 40 --length_required 45 \
      -h  ${outprefix}.fastp.html -j  ${outprefix}.fastp.json -R  ${outprefix}.fsatp.report.txt
else
   for RR in ".1." ".2." "."; do [ -f ${inprefix}${RR}fastq.gz ] && ln -s ${inprefix}${RR}fastq.gz ${outprefix}.clean${RR}fastq.gz ;done
fi
[ $(ls ${outprefix}.clean*fastq.gz 2>/dev/null |wc -l) -gt 0 ] && echo "completed!!" >${snakemake_output[0]}
echo "[## trimmomatic ] date: $(date), parameter: ${sample}, status: finished ..." >>${snakemake_config[status]}
