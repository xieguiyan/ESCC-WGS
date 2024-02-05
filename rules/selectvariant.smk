

rule selectvariant:
   input:
      workp+"/08_Mutect2/{tumor}.part3.completed",
   output:
      completed = workp+"/09_selectvariant/{tumor}.completed",
   threads: 8
   shell:
      """
      outputp=$(echo \"{output.completed}\" |grep -o \"/.*/\")
      if [ ! -d ${{outputp}} ];then mkdir -p ${{outputp}} ;fi
      inputp=$(echo \"{input}\" |grep -o \"/.*/\")
      tumor_vs_normal=$(ls ${{inputp}}|grep -o \"{wildcards.tumor}_vs_[^.]*\"|sort -u)
      sl_inputvcf=\"${{inputp}}/${{tumor_vs_normal}}.filtered0.TN.vcf\"
      if [ $( {java} -jar {GATK} |grep "LeftAlignAndTrimVariants" |wc -l) -gt 0 ]
      then
         time {java} -jar {GATK} LeftAlignAndTrimVariants \\
            -R {reference} \\
            -V ${{inputp}}/${{tumor_vs_normal}}.filtered0.TN.vcf  \\
            -O ${{outputp}}/${{tumor_vs_normal}}.Mt2.biallele.TN.vcf \\
            --split-multi-allelics
         sl_inputvcf=\"${{outputp}}/${{tumor_vs_normal}}.Mt2.biallele.TN.vcf\"
      fi
      time {java} -jar {GATK} SelectVariants \\
         -V ${{sl_inputvcf}} \\
         -O ${{outputp}}/${{tumor_vs_normal}}.Mt2.TN.filter2.vcf   \\
         --exclude-non-variants \\
         --exclude-filtered
      echo "completed!!" >{output.completed}
      echo \"[## filter ] date: $(date), parameter: {wildcards.tumor}, status: finished ...\" >>{status}
      """
