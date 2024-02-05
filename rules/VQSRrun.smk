

rule VQSRrun:
   input:
      Vcf = workp+"/08_haplotypecaller/{sample}.chromall.vcf"
   output:
      recalFile = workp + "/09_vqsrresult/{sample}.snps.recal",
      recalFile2 = workp + "/09_vqsrresult/{sample}.snps.indels.recal",
      tranchesFile = workp + "/09_vqsrresult/{sample}.snps.tranches",
      tranchesFile2 = workp + "/09_vqsrresult/{sample}.snps.indels.tranches",
      rscriptFile = workp + "/09_vqsrresult/{sample}.snps.plots.R",
      rscriptFile2 = workp + "/09_vqsrresult/{sample}.snps.indels.plots.R",
      VQSRvcf = workp + "/09_vqsrresult/{sample}.snps.VQSR.vcf",
      VQSRvcf2 = workp + "/09_vqsrresult/{sample}.snps.indels.VQSR.vcf",
      completed = workp + "/09_vqsrresult/{sample}.completed"
   params:
      Rscript = config["Application"]["Rscript"],
   shell:
      """
      time {java} -jar {GATK}  VariantRecalibrator \\
         -R {reference} -V {input.Vcf} \\
         -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \\
         {bundlehg}/hapmap_3.3.{genomeversion}.vcf.gz \\
         -resource:omini,known=false,training=true,truth=false,prior=12.0 \\
         {bundlehg}/1000G_omni2.5.{genomeversion}.vcf.gz \\
         -resource:1000G,known=false,training=true,truth=false,prior=10.0 {snps1000Gvcf} \\
         -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 {dbsnpvcf} \\
         -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \\
         -mode SNP --dont-run-rscript -O {output.recalFile} \\
         --tranches-file {output.tranchesFile} --rscript-file {output.rscriptFile}

      time {java} -jar {GATK} ApplyVQSR \\
         -R {reference} -V {input.Vcf} -ts-filter-level 99.0 \\
         --tranches-file {output.tranchesFile} --recal-file {output.recalFile} \\
         -mode SNP -O {output.VQSRvcf}

      time {java} -jar {GATK} VariantRecalibrator \\
         -R {reference} --variant {output.VQSRvcf} \\
         -resource:mills,known=true,training=true,truth=true,prior=12.0 {indelsvcf_2} \\
         -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\
         -mode INDEL --dont-run-rscript -O {output.recalFile2} \\
         --tranches-file {output.tranchesFile2} --rscript-file {output.rscriptFile2}

      time {java} -jar {GATK} ApplyVQSR \\
         -R {reference} --variant {output.VQSRvcf} -ts-filter-level 99.0 \\
         --tranches-file {output.tranchesFile2} --recal-file {output.recalFile2} \\
         -mode INDEL -O {output.VQSRvcf2}

      {params.Rscript} {output.rscriptFile}
      {params.Rscript} {output.rscriptFile2}
      echo \"completed!!\" >{output.completed}
      echo \"[## filter ] date: $(date), parameter: {wildcards.sample}, status: finished ...\" >>{status}
      """

