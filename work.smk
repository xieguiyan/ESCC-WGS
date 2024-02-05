##snakemake framework
__Date__ = "Dec 20th, 2023 ~ Jan 5th, 2024"
__Description__ = "ESCC-WGS: WGS analysis for metastatic esophageal squamous cell carcinoma"
__ContactUS__ = "XXXXXXXX@XXX.com"
__copyright__ = "XXXlab"
__license__ = "XXXXXXXX"

##需要修改的参数
workp = "/public/home/caozhou/projects2/PRJNA744411_SRP327447_T"


configfile: workp+"/work.WGS.yaml"
rawdata = config["rawdata"]
status = workp + "/run.status"
config["status"] = status
config["workp2"] = workp
pipeline = config["Pipeline"]


import sys
sys.path.append( pipeline )
from selfpython import DataSimplewrapper
DataSimplewrapper.updatestatus(status)
##因为流程包含了3个大的模块，这里在确定用下面脚本里的哪些规则来实现
resultdict = {"SRA":"03_cleandata", "FASTQ":"03_cleandata", "BAM": "01_rawdata"}
resultnode1 = resultdict[config["Inputfm"]]
NodeandSample = DataSimplewrapper.ObtainNodeandSamplelist(config)
samplelist = NodeandSample["samplelist"]
tumorlist = NodeandSample["tumorlist"]
normallist = NodeandSample["normallist"]
resultnode2, resultnode3 = NodeandSample["resultnode2"], NodeandSample["resultnode3"]
##生成随机字符串，创建临时工作目录，
mkdir_tmpdirectory = DataSimplewrapper.mkdir_tmpdirectory




## 下面开始了流程的rules编写
rule all: 
   input:
      expand( workp+"/" + resultnode2 + "/{sample}.final.completed", sample = samplelist ),
      PPPP = workp+ "/12_finalresult/12_Supplementary.txt"


rule prepare:
   output:
      workp + "/01_rawdata/{sample}.completed"
   params:
      fastqdump = config["Application"]["fastqdump"],
   retries: 3
   log:
      workp + "/01_rawdata/{sample}.log"
   # threads: 8
   script:
      pipeline+"/scripts/prepare.sh"


rule fastqc:
   input:
      workp + "/01_rawdata/{sample}.completed",   
   output:
      workp + "/02_fastqc/{sample}.completed",
   params:
      fastqc = config["Application"]["fastqc"],
      java = config["Application"]["java"],
   resources:
      mem_mb=100,
      # nvidia_gpu=1
   script:
      pipeline+"/scripts/fastqc.sh"


rule cleandata:
   input:
      workp + "/01_rawdata/{sample}.completed",
      completed2 = workp + "/02_fastqc/{sample}.completed",
   output:
      workp+"/03_cleandata/{sample}.completed"
   params:
      fastp = config["Application"]["fastp"],
   # log:
   #   workp+"/03_cleandata/{sample}.log"
   script:
      pipeline+"/scripts/fastprun.sh"


rule BwaRun:
   input:
      workp+"/" + resultnode1 + "/{sample}.completed"
   output:
      BAM = workp+"/04_alignement/{sample}.BWAMEM.bam",
      SortedBam = workp+"/04_alignement/{sample}.BWAMEM.sorted.bam",
   params:
      bwa = config["Application"]["bwa"],
      refer = config["Parameter"]["bwaparams"]["HumanGenomeFasta"],
      samtools = config["Application"]["samtools"],
      marksecondary = DataSimplewrapper.obtainbwapara(config),
      Kvalue = config["bwaKvalue"],
   message:
      "Align lane with bwa mem",
   threads: 8
   script:
      pipeline+"/scripts/bwarun.sh"


rule MarkDuplicate:
   input:
      SortedBam = workp+"/04_alignement/{sample}.BWAMEM.sorted.bam" 
   output:
      SortedMarkdupBam = workp+"/05_markduplicate/{sample}.BWAMEM.sorted.markdup.bam",
      SortedMarkdupBai = workp+"/05_markduplicate/{sample}.BWAMEM.sorted.markdup.bai",
   params:
      java = config["Application"]["java"],
      picard = config["Application"]["picard"],
      samtools = config["Application"]["samtools"],
   script:
      pipeline+"/scripts/markdup.sh"


rule MarkDuplicate_2:
   input:
      SortedBam = workp+"/04_alignement/{sample}.BWAMEM.sorted.bam"
   output:
      SortedMarkdupBam = workp+"/05_markduplicate_2/{sample}.BWAMEM.sorted.markdup.bam",
      SortedMarkdupBai = workp+"/05_markduplicate_2/{sample}.BWAMEM.sorted.markdup.bai",
   params:
      java = config["Application"]["java"],
      samtools = config["Application"]["samtools"],
      sentieon = config["Application"]["sentieon"],
   script:
      pipeline+"/scripts/markdup_2.sh"



#
## 再进行正式的GARK分析之前，处理下跟GATK分析相关的各种参数~~
#
from selfpython import dowithgatkparams
dowithgatkparams.obtainparams(config)
GATK = config["Application"]["GATK"]
java = config["Application"]["java"]
genomeversion = config["genomeversion"]
bundlehg = config["bundlehg"]
indelsvcf_1 = config["indelsvcf_1"]
indelsvcf_2 = config["indelsvcf_2"]
dbsnpvcf = config["dbsnpvcf"]
snps1000Gvcf = config["snps1000Gvcf"]
reference = config["reference"]
afonly = config["afonly"]
small_exac = config["small_exac"]
datatype_gatk = config["datatype_gatk"]
gatkinterval = config["gatkinterval"]
mt2_fir2_need = config["mt2_fir2_need"]
mt2_normal_part1_fomat = config["mt2_normal_part1_fomat"]
mt2_gmr_need = config["mt2_gmr_need"]
mt2_pon_need = config["mt2_pon_need"]
sCNVcallerfinalformat = config["sCNVcallerfinalformat"]
scnv_Allelic_need = config["scnv_Allelic_need"]
mt2_parallerbychrom = config["mt2_parallerbychrom"]
chrom, mt2_chrom, chrom2 = dowithgatkparams.obtainallchromosome(config)
print(config["mt2_fir2_need"], config["mt2_gmr_need"], config["mt2_pon_need"], config["gatkinterval"] )
get_BQSR_bam = dowithgatkparams.get_BQSR_bam
get_normal_realigned_bam = dowithgatkparams.get_normal_realigned_bam


rule alignesentieon: 
   input:
      workp+"/05_markduplicate_2/{sample}.BWAMEM.sorted.markdup.bam"
   output:
      bam = workp+"/06_realignesentieon/{sample}.BWAMEM.sorted.markdup.realigned.bam",
      table = workp+"/06_realignesentieon/{sample}.BWAMEM.sorted.markdup.realigned.table",
   params:
      sentieon = config["Application"]["sentieon"],
      reference = reference,
   threads: 20
   script:
      pipeline+"/scripts/alignesentieon.sh"


rule RealigneTarget:
   input:
      SortedMarkdupBam = workp+"/05_markduplicate/{sample}.BWAMEM.sorted.markdup.bam",
   output:
      ignore = workp+"/06_realigntarget/{sample}.ignore",
   params:
      reference = reference,
      GATK = GATK,
      java = java,
   script:
      pipeline+"/scripts/RealigneTarget.sh"

   
# rule BQSRrun:
#    input:
#       SortedMarkdupBam = workp+"/05_markduplicate/{sample}.BWAMEM.sorted.markdup.bam",
#       ignore = workp+"/06_realigntarget/{sample}.ignore",
#    output:
#       table = workp+"/07_bqsrresult/{sample}.recal_data.table",
#       BqsrBam = workp+"/07_bqsrresult/{sample}.sorted.markdup.BQSR.bam",
#   params:
#       knownsitefile = dowithgatkparams.obtainBQSRknownsite(config)
include: pipeline +"/rules/BQSRrun.smk"


module moduleHaplotypeCaller:
   snakefile: pipeline +"/modules/HaplotypeCaller.module.smk"
   config: config
   skip_validation: True

use rule HaplotypeCallerBychrom from moduleHaplotypeCaller as HCpart1 with:
   input:
      BqsrBam = workp+"/07_bqsrresult/{sample}.sorted.markdup.BQSR.bam"
   output:
      Vcf = workp+"/08_haplotypecaller/{sample}.chromosome/{sample}.{chromno}.part.vcf"
   params:
      reference = reference,
      status = status,
   threads: 24
   resources:
      mem_mb = 5000,
      tmpdir = mkdir_tmpdirectory,
   message:
      "two different method could be chosen for generating VCF..."


use rule HaplotypeCallerAftermerge from  moduleHaplotypeCaller as HCpart2 with:
   input:
      expand(workp+"/08_haplotypecaller/{{sample}}.chromosome/{{sample}}.{chromno}.part.vcf" , chromno = chrom )
   output:
      Vcf = workp+"/08_haplotypecaller/{sample}.chromall.vcf"
   params:
      samtools = config["Application"]["samtools"],
      status = status,


# rule VQSRrun:
#    input:
#       Vcf = workp+"/08_haplotypecaller/{sample}.chromall.vcf"
#    output:
#       recalFile = workp + "/09_vqsrresult/{sample}.snps.recal",
#       recalFile2 = workp + "/09_vqsrresult/{sample}.snps.indels.recal",
#       tranchesFile = workp + "/09_vqsrresult/{sample}.snps.tranches",
#       tranchesFile2 = workp + "/09_vqsrresult/{sample}.snps.indels.tranches",
#       rscriptFile = workp + "/09_vqsrresult/{sample}.snps.plots.R",
#       rscriptFile2 = workp + "/09_vqsrresult/{sample}.snps.indels.plots.R",
#       VQSRvcf = workp + "/09_vqsrresult/{sample}.snps.VQSR.vcf",
#       VQSRvcf2 = workp + "/09_vqsrresult/{sample}.snps.indels.VQSR.vcf",
#       completed = workp + "/09_vqsrresult/{sample}.completed"
#    params:
#       # reference = reference,
#       Rscript=config["Application"]["Rscript"],
include: pipeline+"/rules/VQSRrun.smk"


module modulegCNVcaller:
    snakefile: pipeline+"/modules/gCNVcaller.module.smk"
    config: config
    skip_validation: True


use rule germlinecnvcaller from  modulegCNVcaller with:
   input:
      Vcf = workp + "/09_vqsrresult/{sample}.snps.indels.VQSR.vcf",
      BqsrBam = workp+"/07_bqsrresult/{sample}.sorted.markdup.BQSR.bam",
   output:
      hdf5 = workp + "/10_gcnvcaller/{sample}.BQSR.hdf5",
      completed = workp + "/10_gcnvcaller/{sample}.completed",
   params:
      status = status,


# rule haplotypecallerannotation:
#    input:
#       Vcf = workp + "/09_vqsrresult/{sample}.snps.indels.VQSR.vcf",
#       completed2 = workp + "/10_gcnvcaller/{sample}.completed",
#    output:
#       funcotator = workp+"/11_HaAnnotation/{sample}.VQSR.funcotator.vcf",
#       completed = workp+"/11_HaAnnotation/{sample}.final.completed",
include: pipeline+"/rules/harannotation.smk"


module modulemutec2:
   snakefile: pipeline +"/modules/Mutect2.module.smk"
   config: config
   skip_validation: True


use rule Mutect2part1ignore from modulemutec2 with:
   input:
     BqsrBam = workp+"/07_bqsrresult/{normal}.sorted.markdup.BQSR.bam",
   output:
     ignore = workp+"/08_Mutect2/{normal}.Mt2.normal.ignore",


use rule Mutect2part1Bychrom from modulemutec2 with:
   input:
      BqsrBam = workp+"/07_bqsrresult/{normal}.sorted.markdup.BQSR.bam",
   output:
      vcfpart = workp+"/08_Mutect2/{normal}.chromosome.vcf/{normal}.{chromno}.Mt2.normal.part.vcf",
   resources:
      tmpdir = mkdir_tmpdirectory,
      mem_mb = 5000,
   threads: 24
   params:
      status = status,


use rule Mutect2part1 from modulemutec2 with:
    input:
       expand(workp+"/08_Mutect2/{normal}.chromosome.vcf/{normal}.{chromno}.Mt2.normal.part.vcf", normal = normallist, chromno = mt2_chrom)
    output:
       expand(workp+"/08_Mutect2/{normal}.Mt2.normal.vcf", normal = normallist )
       # pon = workp+"/08_Mutect2/wholeNoraml.Mt2.pon.vcf",
    threads: 8
    params:
       status = status,


use rule Mutect2part2 from modulemutec2 with:
   input:
      normal = lambda wildcards: workp+"/07_bqsrresult/"+ config["tumor_vs_normal"][wildcards.tumor] +".sorted.markdup.BQSR.bam",
      tumor = workp+"/07_bqsrresult/{tumor}.sorted.markdup.BQSR.bam",
   output:
      completed = workp+"/08_Mutect2/{tumor}.part2.1.completed"
   threads: 8
   params:
      status = status,


use rule Mutect2part2bychrom from modulemutec2 with:
   input:
      Mtnormal = lambda wildcards: workp+"/08_Mutect2/" + config["tumor_vs_normal"][wildcards.tumor] + ".Mt2.normal." + mt2_normal_part1_fomat,
      tumor = workp+"/07_bqsrresult/{tumor}.sorted.markdup.BQSR.bam",
      completed = workp+"/08_Mutect2/{tumor}.part2.1.completed",
   output:
      completed = workp+"/08_Mutect2/{tumor}.Mutect2.chromosome/{tumor}.{chromno}.part2.2.completed",
   resources:
      tmpdir = mkdir_tmpdirectory,
      mem_mb = 5000,
   threads: 24
   params:
      status = status,


use rule Mutect2part3 from modulemutec2 with:
   input:
      expand(workp+"/08_Mutect2/{{tumor}}.Mutect2.chromosome/{{tumor}}.{chromno}.part2.2.completed", chromno = mt2_chrom )
   output:
      completed = workp+"/08_Mutect2/{tumor}.part3.completed"
   threads: 8
   params:
      status = status,

   
# rule selectvariant:
#    input:
#       workp+"/08_Mutect2/{tumor}.part3.completed",
#    output:
#       completed = workp+"/09_selectvariant/{tumor}.completed",
#    threads: 8
include: pipeline + "/rules/selectvariant.smk"


module modulesCNVcaller:
   snakefile: pipeline+"/modules/sCNVcaller.module.smk"
   config: config
   skip_validation: True


use rule gatkCNVcallerignore from modulesCNVcaller with:
   input:
      normal = lambda wildcards: get_BQSR_bam(wildcards, config, workp),
      tumor = workp+"/07_bqsrresult/{tumor}.sorted.markdup.BQSR.bam",
   output:
      ignore = workp+"/10_cnvcaller/{tumor}.part2.2.ignore"


use rule gatkCNVcallerpairpart1 from modulesCNVcaller with:
   input:
      expand(workp+"/07_bqsrresult/{sample}.sorted.markdup.BQSR.bam", sample = samplelist )
      # normal = lambda wildcards: get_BQSR_bam(wildcards, config, workp),
      # tumor = workp+"/07_bqsrresult/{tumor}.sorted.markdup.BQSR.bam",
   output:
      expand( workp+"/10_cnvcaller/{tumor}.part1.completed", tumor = tumorlist )
      # completed = workp+"/10_cnvcaller/{tumor}.part1.completed",
   resources:
      mem_mb=10000 # 就在这里报错了内存溢出，堆栈错误的情况，所以就特意分配了内存空间
   params:
      status = status,


use rule gatkCNVcallerpairpart2bychrom from modulesCNVcaller with:
   input:
      tumor = workp+"/07_bqsrresult/{tumor}.sorted.markdup.BQSR.bam",
      completed = workp+"/10_cnvcaller/{tumor}.part1.completed",
   output:
      completed = workp+"/10_cnvcaller/{tumor}.chromosome/{tumor}.{chromno}.part2.1.completed",
   resources:
      mem_mb = 10000,
      tmpdir = mkdir_tmpdirectory,
   params:
      normal = lambda wildcards: config["tumor_vs_normal"][wildcards.tumor],
      status = status,


use rule gatkCNVcallerpairpart2merge from modulesCNVcaller with:
   input:
      expand( workp+"/10_cnvcaller/{{tumor}}.chromosome/{{tumor}}.{chromno}.part2.1.completed", chromno = chrom2 )
   output:
      completed = workp+"/10_cnvcaller/{tumor}.part2.2.completed",
   params:
      tumor_vs_normal = lambda wildcards: wildcards.tumor + "_vs_"+ config["tumor_vs_normal"][wildcards.tumor],
      GATKsCNVchrom = config["Parameter"]["GATKparams"]["sCNVchrom"],
   resources:
      mem_mb = 20000,
      tmpdir = mkdir_tmpdirectory,



# rule mutect2annote:
#    input:
#       workp+"/08_Mutect2/{tumor}.part3.completed",
#       completed = workp+"/09_selectvariant/{tumor}.completed",
#       completed2 = workp+"/10_cnvcaller/{tumor}.part2.2." + sCNVcallerfinalformat,
#    output:
#       completed = workp+"/11_mutect2annotation/{tumor}.completed",
include: pipeline + "/rules/mutect2annote.smk"
   

module moduleECGEA:
   snakefile: pipeline+"/modules/ECGEA.module.smk"
   config: config
   skip_validation: True


use rule TNhaplotyper  from  moduleECGEA with:
   input:
      normal = lambda wildcards: get_normal_realigned_bam(wildcards, config, workp),
      tumor = workp+"/06_realignesentieon/{tumor}.BWAMEM.sorted.markdup.realigned.bam",
      tumortb = workp+"/06_realignesentieon/{tumor}.BWAMEM.sorted.markdup.realigned.table",
   output:
      completed = workp+"/07_TNhaplotyper/{tumor}.completed",
   params:
      sentieon = config["Application"]["sentieon"],
      bedtools = config["Application"]["bedtools"],
      reference = reference,
   threads: 8
   params:
      status = status,


use rule cnvkitrun from  moduleECGEA with:
   input:
      normal = get_normal_realigned_bam,
      tumor = workp+"/06_realignesentieon/{tumor}.BWAMEM.sorted.markdup.realigned.bam",
      completed = workp+"/07_TNhaplotyper/{tumor}.completed"
   output:
      completed = workp + "/08_cnvkitresult/{tumor}.completed"
   params:
      cnvkit = config["Application"]["cnvkit"],
      python = config["Application"]["python"],
      reference = reference,
   params:
      status = status,


# rule snpeffrun:
#    input:
#       completed = workp+"/07_TNhaplotyper/{tumor}.completed",
#       completed2 = workp + "/08_cnvkitresult/{tumor}.completed",
#    output:
#       completed = workp + "/09_ECannotation/{tumor}.completed",
#    params:
#       snpEff = config["Application"]["snpEff"],
#       annote = lambda wildcards: workp + "/09_ECannotation/%s_vs%s.TN.filter2.annotate.vcf"%(wildcards.tumor, config["tumor_vs_normal"][wildcards.tumor])
include: pipeline + "/rules/snpeffrun.smk"


rule AfterTumorVsNormal:
   input:
      expand( workp+"/"+resultnode3+"/{tumor}.completed", tumor = tumorlist)
   output:
      expand( workp+"/" + resultnode3 + "/{sample}.final.completed", sample = samplelist )
   script:
      pipeline+"/scripts/afterTN.sh"


# rule Suplementary:
#    input:
#       expand( workp+"/" + resultnode2 + "/{sample}.final.completed", sample = samplelist )
#    output:
#       PPPP = workp+ "/12_finalresult/12_Supplementary.txt",
#    params:
#       Supple = config["Files"]["Supplementary"],
#       MainModule = config["MainModule"]
include: pipeline+"/rules/supplementary.smk"
