
configfile: "null.yaml"

GATK = config["Application"]["GATK"]
java = config["Application"]["java"]
genomeversion = config["genomeversion"]
bundlehg = config["bundlehg"]
indelsvcf_1 = config["indelsvcf_1"]
indelsvcf_2 = config["indelsvcf_2"]
dbsnpvcf = config["dbsnpvcf"]
reference = config["reference"]



workp = ""
chrom = ["chr1"]



rule TNhaplotyper:
   input:
      normal = "null",
      tumor = workp+"/06_realignesentieon/{tumor}.BWAMEM.sorted.markdup.realigned.bam",
      tumortb = workp+"/06_realignesentieon/{tumor}.BWAMEM.sorted.markdup.realigned.table",
   output:
      completed = workp+"/07_TNhaplotyper/{tumor}.completed",
   params:
      sentieon = config["Application"]["sentieon"],
      bedtools = config["Application"]["bedtools"],
   threads: 8
   params:
      status = "null",
   shell:
      """
      normaltb=$(echo \"{input.normal}\"|sed \"s/bam$/table/\")
      normal=$(echo \"{input.normal}\"|grep \"[^/]*$\" |cut -d \".\" -f1)
      tumor=$(echo \"{input.tumor}\"|grep \"[^/]*$\" |cut -d \".\" -f1)
      outputp=$(echo \"{output.completed}\" | grep -o \"/.*/\")
      if [ ! -d ${{outputp}} ]; then mkdir -p ${{outputp}} ;fi
      prefix=${{outputp}}/${{tumor}}_vs_${{nromal}}
      sentieondata=$(echo \"{params.sentieon}\"|sed \"s/bin\/.*$//\")/data
      time {params.sentieon} driver -t {threads} \\
          -i {input.tumor} -q {input.tumortb} -i {input.normal} -q ${{normaltb}} \\
          -r {reference} --algo TNhaplotyper \\
          --tumor_sample {wildcards.tumor} --normal_sample ${{normal}} \\
          --dbsnp {dbsnpvcf} --cosmic ${{sentieondata}}/CosmicCodingMuts.normal.vcf \\
          ${{prefix}}.TN.vcf
      
      #{params.bedtools} subtract -header -a `# ${{tumor}}_vs_${{normal}}_TN.vcf` {{output.TNVCF}} \\
      #    -b ${{black_region}} | awk -F \'\t\' \'($7==\"PASS\" || $7==\"clustered_events\") || \\
      #       $1~\"^#\"\' | uniq > ${{prefix}}.TN.PASS.vcf  #{{output.TNPASS}}

      cat ${{prefix}}.TN.vcf | awk -F \'\t\' \'($7==\"PASS\" || $7==\"clustered_events\") || \
         $1~\"^#\"\' | uniq > ${{prefix}}.TN.PASS.vcf  #{{output.TNPASS}}
      echo \"completed!!\" >{output.completed}
      echo \"[## TNhaplotyper ] date: $(date), parameter: {wildcards.tumor}, status: finished ...\" >>{params.status}
      """


rule cnvkitrun:
   input:
      normal = "null",
      tumor = workp+"/06_realignesentieon/{tumor}.BWAMEM.sorted.markdup.realigned.bam",
      completed = workp+"/07_TNhaplotyper/{tumor}.completed"
   output:
      completed = workp + "/08_cnvkitresult/{tumor}.completed"
   params:
      cnvkit = config["Application"]["cnvkit"],
      python = config["Application"]["python"],
   params:
      status = "null",
   shell:
      """
      cnvkitp=$(echo \"{params.cnvkit}\" |grep -o \"/.*/\")
      outputp=$(echo \"{output}\" |grep -o \"/.*/\" )
      if [ ! -d ${{outputp}} ]; mkdir -p ${{outputp}} ;fi
      tumor=$(echo \"{input.tumor}\"|grep -o \"[^/]*$\" |cut -d \".\" -f1)
      normal=$(echo \"{input.normal}\"|grep -o \"[^/]*$\" |cut -d \".\" -f1)
      tumor_vs_normal=${{tumor}}_vs_${{normal}}
      {params.python} {params.cnvkit} batch {input.tumor} --normal {input.normal} \\
         `# --targets my_baits.bed` --method wgs --fasta {reference} \\
         `# --access  data/access-5kb-mappable.hg19.bed` \\
         --access ${{cnvkitp}}/../data/access-10kb.{genomeversion}.bed \\
         `# --antitergets panel.antitarget.bed` \\
         --annotate ${{cnvkitp}}/../data/refFlat_{genomeversion}.txt \\
         --output-reference ${{outputp}}/${{tumor_vs_normal}}_reference.cnn \\
         --output-dir ${{outputp}}
      {params.python} {params.cnvkit} batch {input.tumor} \\
         -r ${{outputp}}/${{tumor_vs_normal}}_reference.cnn \\
         -p 0 --scatter --diagram -d ${{outputp}}/
      echo \"completed!!\" >{output.completed}
      echo \"[## sCNVcaller ] date: $(date), parameter: {wildcards.tumor}, status: finished ...\" >>{params.status}
      """

