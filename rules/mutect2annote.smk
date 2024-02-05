


rule mutect2annote:
   input:
      workp+"/08_Mutect2/{tumor}.part3.completed",
      completed = workp+"/09_selectvariant/{tumor}.completed",
      # completed2 = workp+"/10_cnvcaller/{tumor}.part2.2.completed",
      completed2 = workp + "/10_cnvcaller/{tumor}.part2.2." + sCNVcallerfinalformat,
   output:
      completed = workp+"/11_mutect2annotation/{tumor}.completed",
   shell:
      """
      inputp=$(echo \"{input.completed}\"|grep -o \"/.*/\")
      tumor_vs_normal=$(ls ${{inputp}}| grep -o \"{wildcards.tumor}_vs_[^.]*\"|sort -u)
      filter_vcf=${{inputp}}${{tumor_vs_normal}}.Mt2.TN.filter2.vcf
      raw_vcf=$(echo \"{input[0]}\" |grep -o \"/.*/\" )${{tumor_vs_normal}}.Mt2.TN.vcf
      outputp=$(echo \"{output.completed}\" |grep -o \"/.*/\" )
      funcotatorpath={bundlehg}/../funcotator
      for vcfone in \"${{filter_vcf}}\" \"${{raw_vcf}}\"
      do
         time {java} -jar {GATK} Funcotator \\
            -R {reference} -V ${{vcfone}} \\
            -O ${{outputp}}/$(echo \"${{vcfone}}\" |sed \"s/\/.*\///\" |sed \"s/.vcf$//\" ).funcotator.vcf \\
            --data-sources-path \\
            ${{funcotatorpath}}/funcotator_dataSources.v1.8.{genomeversion}.20230908s/ \\
            --ref-version {genomeversion} --output-file-format VCF
      done
      echo \"completed!!\" >{output.completed}
      echo \"[## annotate ] date: $(date), parameter: {wildcards.tumor}, status: finished ...\" >>{status}
      """

