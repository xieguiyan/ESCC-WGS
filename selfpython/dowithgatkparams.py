#!/usr/bin/python
import os
import re


def TrueorFalse(para):
  if isinstance(para, bool) and para: return "true"
  elif isinstance(para,str) and para[0].lower()== "t": return "true"
  else: return "false"


def obtainparams(config):
   GATK = config["Application"]["GATK"]
   java = config["Application"]["java"]
   genomeversion = config["Parameter"]["GATKparams"]["genomeversion"]
   config["genomeversion"] = genomeversion
   bundle = config["Parameter"]["GATKparams"]["bundle"]
   bundlehg = bundle + "/" + genomeversion
   config["bundlehg"] = bundlehg
   indelsvcf_1 = bundlehg + "/1000G_phase1.indels.%s.vcf.gz"%genomeversion
   config["indelsvcf_1"] = indelsvcf_1
   indelsvcf_2 = bundlehg + "/Mills_and_1000G_gold_standard.indels.%s.vcf.gz"%genomeversion
   config["indelsvcf_2"] = indelsvcf_2
   dbsnpvcf = bundlehg + "/dbsnp_" + str( config["Parameter"]["GATKparams"]["dbsnp_v"] ) + f".{genomeversion}.vcf.gz"
   config["dbsnpvcf"] = dbsnpvcf
   snps1000Gvcf = bundlehg + "/1000G_phase1.snps.high_confidence.%s.vcf.gz"%genomeversion
   config["snps1000Gvcf"] = snps1000Gvcf
   if "hg38" in genomeversion:
      # reference = bundlehg + f"/Homo_sapiens_assembly{genomeversion[2:]}.fasta"
      reference = bundlehg +f"/bwa_index_BuildByMyself2/hg38.fa"
   elif "hg19" in genomeversion:
      reference = bundlehg + "/ucsc.hg19.fasta"
   else:
      print(" ERROR: judge genome version failed!!")
      exit(-1)
   config["reference"] = reference
   afonly = f"{bundlehg}/../Mutect2/af-only-gnomad.{genomeversion}.vcf.gz"
   config["afonly"] = afonly
   small_exac = "%s/../Mutect2/GetPileupSummaries/small_exac_common_3.%s.vcf.gz"%(bundlehg, genomeversion )
   config["small_exac"] = small_exac
   
   datatype_gatk = config["Parameter"]["GATKparams"]["datatype"]
   config["datatype_gatk"] = datatype_gatk
   if config["Parameter"]["GATKparams"]["interval"] == "NULL.bed":
      gatkinterval = bundlehg + "/wgs_calling_regions.%s.interval_list"%genomeversion
   else:
      gatkinterval = config["Parameter"]["GATKparams"]["interval"]
   config["gatkinterval"] = gatkinterval
   mt2_fir2_need = TrueorFalse( config["Parameter"]["GATKparams"]["Mutect2params"]["fir2_need"] )
   config["mt2_fir2_need"] = mt2_fir2_need
   mt2_normal_part1_fomat = "vcf"
   if mt2_fir2_need == "false" : mt2_normal_part1_fomat = "ignore"
   config["mt2_normal_part1_fomat"] = mt2_normal_part1_fomat
   mt2_gmr_need = TrueorFalse( config["Parameter"]["GATKparams"]["Mutect2params"]["gmr_need"] )
   config["mt2_gmr_need"] = mt2_gmr_need
   mt2_pon_need = TrueorFalse( config["Parameter"]["GATKparams"]["Mutect2params"]["pon_need"] )
   config["mt2_pon_need"] = mt2_pon_need
   mt2_parallerbychrom = config["Parameter"]["GATKparams"]["Mutect2params"]["parallerbychrom"]
   config["mt2_parallerbychrom"] = mt2_parallerbychrom
   sCNVcallerrun = TrueorFalse( config["Parameter"]["GATKparams"]["sCNVcaller"]["run_need"] )
   config[ "sCNVcallerrun" ] = sCNVcallerrun
   if sCNVcallerrun:
      sCNVcallerfinalformat = "ignore"
   else:
      sCNVcallerfinalformat = "completed"
   config["sCNVcallerfinalformat"] = sCNVcallerfinalformat
   config["Parameter"]["GATKparams"]["sCNVchrom"] = config["Parameter"]["GATKparams"]["sCNVcaller"]["bychrom"]
   config["scnv_Allelic_need"] =  TrueorFalse( config["Parameter"]["GATKparams"]["sCNVcaller"]["Allelic_need"] )
  
   

def obtainBQSRknownsite(config):
   vcffile=""
   for vcfone in config["Parameter"]["GATKparams"]["BQSRparams"]["known_sites"]:
      vcffile += " " + config[vcfone]
   return vcffile   

def obtainallchromosome(config):
   chrom = []
   for number in range(1, 23): chrom.append("chr%d"%number)
   chrom += ["chrM","chrX", "chrY"]
   mt2_chrom = ["chromall"]
   if ( config["mt2_parallerbychrom"] != "chromall" ):
      mt2_chrom = chrom
   chrom2 = chrom if config["Parameter"]["GATKparams"]["sCNVchrom"] == "bychrom" else [ "chromall" ]
   config["chrom"] = chrom
   return [chrom, mt2_chrom, chrom2 ]

def get_BQSR_bam(wildcards, config, workp):
   return workp+"/07_bqsrresult/"+ config["tumor_vs_normal"][wildcards.tumor] +".sorted.markdup.BQSR.bam"

def get_normal_realigned_bam(wildcards, config, workp):
   prefix =  workp+"/06_realignesentieon/"+ config["tumor_vs_normal"][wildcards.tumor]
   return prefix + ".BWAMEM.sorted.markdup.realigned.bam"
