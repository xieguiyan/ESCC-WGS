#!/bin/bash
sample="${snakemake_wildcards[sample]}"
inputp=$(echo "${snakemake_input[0]}" |grep -o "/.*/")
outputp=$(echo "${snakemake_output[0]}" |sed "s/[^/]*$//" )
if [ ! -d ${outputp} ] ;then mkdir -p ${outputp} ;fi
if [[ $(ls ${inputp} |grep "${sample}.*\.fastq\.gz" |wc -l) -gt 0 ]]
then
   ${snakemake_params[fastqc]} -o ${outputp} -j \
       ${snakemake_params[java]} ${inputp}/${sample}*.fastq.gz
   QT="PASS" #base quality
   QT_criterion1="Per sequence quality scores"
   QT_criterion1="Adapter Content"
   for htmlfile in $(ls ${outputp}/${sample}*_fastqc.html )
   do
      before_that="alt\=\"\[PASS\]\".{5,6}href\=\"\#M[0-9]{1,2}\">"
      grep -E -o "${before_that}${QT_criterion1}" $htmlfile
      if [ $(grep -E -o "${before_that}${QT_criterion1}" $htmlfile |wc -l) -lt 1 ] || \
         [ $(grep -E -o "${before_that}${QT_criterion2}" $htmlfile |wc -l) -lt 1 ]
      then
         QT="FAILED" ;continue
      fi
   done
   echo "extract result from fastqc..." >${outputp}/${sample}.fastqc.${QT}
fi
[ $(ls ${outputp}/${sample}*_fastqc.html 2>/dev/null|wc -l) -gt 0 ] && echo "completed!!" >${snakemake_output[0]}
echo "[## fastqc ] date: $(date), parameter: ${sample}, status: finished ..." >>${snakemake_config[status]}
