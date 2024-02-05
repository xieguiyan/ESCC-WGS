

rule haplotypecallerannotation:
   input:
      Vcf = workp + "/09_vqsrresult/{sample}.snps.indels.VQSR.vcf",
      completed2 = workp + "/10_gcnvcaller/{sample}.completed",
   output:
      funcotator = workp+"/11_HaAnnotation/{sample}.VQSR.funcotator.vcf",
      completed = workp+"/11_HaAnnotation/{sample}.final.completed",
   shell:
      """
      # time {java} -jar {GATK} CNNScoreVariants \\
      #   -V {{input.Vcf}} -R {reference} -O ${{output.CNNvcf}}

      #time {java} -jar {GATK} FilterVariantTranches \\
      #   -V ${{output.CNNvcf}} \\
      #   --resource {indelsvcf_2} \\
      #   --resource {dbsnpvcf} \\
      #   --resource {snps1000Gvcf} \\
      #   --resource {bundlehg}/1000G_omni2.5.{genomeversion}.vcf \\
      #   --resource {bundlehg}/hapmap_3.3.{genomeversion}.vcf.gz \\
      #   --info-key CNN_1D --snp-tranche 99.0 \\
      #   -O ${{output.filter2}}
      
      funcotatorpath={bundlehg}/../funcotator
      time {java} -jar {GATK} Funcotator \\
         -R {reference} -V {input.Vcf} \\
         -O {output.funcotator} \\
         --data-sources-path \\
         ${{funcotatorpath}}/funcotator_dataSources.v1.8.{genomeversion}.20230908g/ \\
         --ref-version {genomeversion} --output-file-format VCF
      echo \"completed!!\" >{output.completed}
      echo \"[## annotate ] date: $(date), parameter: {wildcards.sample}, status: finished ...\" >>{status}
      """
