
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


rule germlinecnvcaller:
   input:
      Vcf = workp + "/09_vqsrresult/{sample}.snps.indels.VQSR.vcf",
      BqsrBam = workp+"/07_bqsrresult/{sample}.sorted.markdup.BQSR.bam",
   output:
      hdf5 = workp + "/10_gcnvcaller/{sample}.BQSR.hdf5",
      completed = workp + "/10_gcnvcaller/{sample}.completed",
   # log:
   #    # workp + "/10_gcnvcaller/{sample}.log",
   params:
      status = "null",
   shell:
       """
       echo \"待测试的源码，可以参照， gCNVcaller.module.smk.bk\"
       echo \"null, Actually this file shoule be a binary file, skip...\" > {output.hdf5}
       echo \"completed!!\" >{output.completed}
       echo \"[## gCNVcaller ] date: $(date), parameter: {wildcards.sample}, status: finished ...\" >>{params.status}
       """

