#!/bin/bash

sstatus="${snakemake_config[status]}"
sample="${snakemake_wildscards[sample]}"
sentieon="${snakemake_params[sentieon]}"
echo "[## realigned ] date: $(date), parameter: ${sample}, status: running ..." >>${sstatus}
indel_1kg_para="-k ${snakemake_config[indelsvcf_1]}"
if [ ! -f ${snakemake_config[indelsvcf_1]} ]
then 
   indel_1kg_para=""
fi
time ${sentieon} driver -t ${snakemake[threads]} \
   -r ${snakemake_params[reference]} \
   -i ${snakemake_input[0]} --algo Realigner \
   -k `# ${indel_mills}}` ${snakemake_config[indelsvcf_2]} \
   ${indel_1kg_para} `# -k  ${indel_1kg}` \
   ${snakemake_output[bam]}
time ${sentieon} driver -t ${snakemake[threads]} \
   -r ${snakemake_params[reference]} \
   -i ${snakemake_output[bam]} --algo QualCal \
   ${indel_1kg_para} `# -k ${indel_1kg}` \
   -k `# ${indel_mills}` ${snakemake_config[indelsvcf_2]} \
   -k `# ${dbsnp}` ${snakemake_config[dbsnpvcf]} \
   ${snakemake_output[table]}
echo "[## realigned ] date: $(date), parameter: ${sample}, status: finished ..." >>${sstatus}

