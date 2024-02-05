
configfile: "null.yaml"

GATK = config["Application"]["GATK"]
java = config["Application"]["java"]
genomeversion = config["genomeversion"]
bundlehg = config["bundlehg"]
indelsvcf_1 = config["indelsvcf_1"]
indelsvcf_2 = config["indelsvcf_2"]
dbsnpvcf = config["dbsnpvcf"]
reference = config["reference"]
workp = config["workp2"]
chrom = config["chrom"]

rule HaplotypeCallerBychrom:
   input:
      BqsrBam = workp+"/07_bqsrresult/{sample}.sorted.markdup.BQSR.bam"
   output:
      Vcf = workp+"/08_haplotypecaller/{sample}.chromosome/{sample}.{chromno}.part.vcf"
   params:
      reference = reference,
      status = "run.status",
   threads: 16
   resources:
      mem_mb = 5000,
      tmpdir = "null",
   message:
      "two different method could be chosen for generating VCF..."
   shell:
      """
      UsedMethod=1 #1,2
      gVcf=$(echo {output.Vcf} |grep -o \"/.*\.\")g.vcf
      bamout=$(echo {output.Vcf} |grep -o \"/.*\.\")bamout.bam
      if [ ${{UsedMethod}} -eq 1 ]
      then
         time {java} -Djava.io.tmpdir={resources.tmpdir} -jar {GATK} HaplotypeCaller \\
            -R {reference} -I {input.BqsrBam} --emit-ref-confidence GVCF \\
            -L {wildcards.chromno} --bam-output ${{bamout}} -O ${{gVcf}}
         time {java} -jar {GATK} GenotypeGVCFs \\
            -R {reference} -V ${{gVcf}} -O {output.Vcf}

      elif [ ${{UsedMethod}} -eq 2 ]
      then
         time {java} -jar {GATK} HaplotypeCaller \\
            -R {reference} -I {input.BqsrBam} `# -D {dbsnpvcf}` \\
            `# -stand_call_conf 50 -A QualByDepth` \\
            -L {wildcards.chromno} `# -A RMSMappingQuality` `# -A MappingQualityRankSumTest` \\
            `# -A ReadPosRankSumTest` `# -A FisherStrand` `# -A StrandOddsRatio` \\
            `# -A Coverage` -O {output.Vcf}
      fi
      """


rule HaplotypeCallerAftermerge:
   input:
      expand(workp+"/08_haplotypecaller/{{sample}}.chromosome/{{sample}}.{chromno}.part.vcf" , chromno = chrom )
   output:
      Vcf = workp+"/08_haplotypecaller/{sample}.chromall.vcf"
   params:
      samtools = config["Application"]["samtools"],
      status = "run.status",
   shell:
      """
      mergeline=\"\"
      for vcfone in {input};do mergeline=\"${{mergeline}} -I ${{vcfone}}\" ;done
      {java} -jar {GATK} MergeVcfs ${{mergeline}} -O {output.Vcf}
      bamout=$(echo {output.Vcf} |grep -o \"/.*\.\")bamout.bam
      {params.samtools} merge -f ${{bamout}} $(echo \"{input}\" |sed \"s/part\.vcf/part\.bamout.bam/g\")
      # bgzip -f {output.Vcf}
      # tabix -p vcf {output.Vcf}
      echo \"[## haplotypecaller ] date: $(date), parameter: {wildcards.sample}, status: finished ...\" >>{params.status}
      """
