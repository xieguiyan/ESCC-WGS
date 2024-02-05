#!/bin/bash
sample="${snakemake_wildcards[sample]}"
rawdata="${snakemake_config[rawdata]}"
wkrawdata=$( echo "${snakemake_output[0]}" |grep -o "/.*/" )
if [ ! -d ${wkrawdata} ] ;then mkdir -p ${wkrawdata} ;fi
for datain in  "${rawdata}/${sample}" "${rawdata}"
do
   if [[ ${snakemake_config[Inputfm]} == "BAM" ]] && [ -f ${datain}/${sample}.bam ] && \
      [ ! -f ${wkrawdata}/${sample}.bam ]
   then
      ln -s ${datain}/${sample}.bam ${wkrawdata}/ 1>${snakemake_log[0]}
   elif [ $(ls ${datain}/ |grep "${sample}.*\.fastq\.gz$" |wc -l ) -gt 0 ] && \
      [ $(ls  ${wkrawdata}/ | grep "${sample}.*fastq\.gz$" |wc -l ) -lt 1 ]
   then
      ln -s ${datain}/${sample}*.fastq.gz ${wkrawdata}/ 1>${snakemake_log[0]}
   elif [ -f ${datain}/${sample}.sra ] && [ ! -f  ${wkrawdata}/${sample}.sra ]
   then
      ln -s ${datain}/${sample}.sra ${wkrawdata}/
   fi
done
if [ -f ${wkrawdata}/${sample}.sra ] && \
   [ $(ls  ${wkrawdata}/ |grep "${sample}.*fastq.gz" |wc -l ) -lt 1 ] 
then 
   ${snakemake_params[fastqdump]} --split-3 --gzip\
      -O ${wkrawdata}/ ${wkrawdata}/${sample}.sra 1>${snakemake_log[0]}
fi

if [ $(ls ${wkrawdata}/${sample}*fastq.gz 2>/dev/null |wc -l) -gt 0 ]
then 
   echo "completed!!" >${snakemake_output[0]}
   for fastq in  $(ls ${wkrawdata}/${sample}*fastq.gz)
   do
      # Name1=$(echo "$fastq" |grep -E -o "[_rR]+[12].fastq.gz$")
      Name1=$(echo "$fastq" |grep -E -o "[\._rR][12].fastq.gz$")
      Name2=${Name1:1}
      [ -f ${wkrawdata}/${sample}.${Name2} ] && continue
      mv $fastq ${wkrawdata}/${sample}.${Name2}
   done
fi
echo "[## prepare ] date: $(date), parameter: ${sample}, status: finished ..." >>${snakemake_config[status]}
